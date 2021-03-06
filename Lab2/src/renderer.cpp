#include "renderer.h"

#include "application.h"
#include "camera.h"
#include "shader.h"
#include "mesh.h"
#include "texture.h"
#include "fbo.h"
#include "prefab.h"
#include "material.h"
#include "utils.h"
#include "scene.h"
#include "extra/hdre.h"

#include <iostream>
#include <algorithm>
#include <vector>  
#include <cmath>


using namespace GTR;


GTR::Renderer::Renderer()
{
	max_lights = 10;

	// FBOs
	fbo = NULL;
	gbuffers_fbo = NULL;
	illumination_fbo = NULL;
	shadowmap = NULL;
	ssao_fbo = NULL;
	screen_texture = NULL;
	screen_fbo = NULL;
	width_shadowmap = 2048;
	height_shadowmap = 2048;

	// Debug parameters
	show_shadowmap = false;
	debug_shadowmap = 7;
	debug_texture = eTextureType::COMPLETE;
	show_buffers = false;
	pbr = true;

	pipeline = DEFERRED;
	typeOfRender = eRenderPipeline::MULTIPASS;

	// Color correction
	gamma = true;
	tonemapping = true;
	lumwhite2 = 0.8;
	averagelum = 1.6;
	scale = 1.0;

	// Ambient occlusion
	random_points_sph = generateSpherePoints(64, 1, false);
	random_points_hemi = generateSpherePoints(64, 1, true);
	show_ssao = false;
	SSAOType = SSAO_plus;

	probe.pos.set(82, 92, -2 29);

	// meshes
	sphere = Mesh::Get("data/meshes/sphere.obj", false, false);
	quad = Mesh::getQuad();
}

// --- Rendercalls manager functions ---

// Generate the rendercalls vector by iterating through the entities vector
void GTR::Renderer::createRenderCalls(GTR::Scene* scene, Camera* camera)
{
	render_calls.clear();

	// Iterate the entities vector to save each node
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];

		// Save only the visible nodes
		if (!ent->visible)
			continue;

		// If prefab iterate the nodes
		if (ent->entity_type == PREFAB)
		{
			PrefabEntity* pent = (GTR::PrefabEntity*)ent;
			// First take the root node
			if (pent->prefab) {
				GTR::Node* curr_node = &pent->prefab->root;
				Matrix44 node_model = curr_node->getGlobalMatrix(true) * ent->model;
				// Pass the global matrix for this node and the root one to compute the global matrix for the rest of children
				addRenderCall_node(scene, camera, curr_node, node_model, ent->model);
			}
		}
	}
}

// Recursive function to add a rendercall node with its children
void GTR::Renderer::addRenderCall_node(GTR::Scene* scene, Camera* camera, Node* node, Matrix44 curr_model, Matrix44 root_model) {
	RenderCall rc;
	// If the node doesn't have mesh or material do not add it
	if (node->material && node->mesh) {
		Vector3 nodepos = curr_model.getTranslation();
		rc.mesh = node->mesh;
		rc.material = node->material;
		rc.model = curr_model;
		rc.distance_to_camera = nodepos.distance(scene->main_camera.eye);
		rc.world_bounding = transformBoundingBox(curr_model, node->mesh->box);
		// If the material is opaque add a distance factor to sort it at the end of the vector
		if (rc.material->alpha_mode == GTR::eAlphaMode::BLEND)
		{
			int dist_factor = 1000000;
			rc.distance_to_camera += dist_factor;
		}
		render_calls.push_back(rc);
	}


	// Add also all the childrens of this node
	for (int j = 0; j < node->children.size(); ++j) {
		GTR::Node* curr_node = node->children[j];
		// Compute global matrix
		Matrix44 node_model = node->getGlobalMatrix(true) * root_model;
		addRenderCall_node(scene, camera, curr_node, node_model, root_model);
	}
}

// Sort rendercalls by distance
void GTR::Renderer::sortRenderCalls() {
	std::sort(render_calls.begin(), render_calls.end(), compare_distances);
}

// --- Shadowmap functions ---

// to assign row, column number and size of the grid in an efficient way specially for shadowmap
Vector4 GTR::Renderer::assignMapPiece(int width, int height, int index, int num_elements) {
	// get number of cols and rows we need for the number of elements we need to store
	int num_cols = (int)ceil(sqrt(num_elements));
	int num_rows;
	// if last row remains empty delete it
	if (num_elements <= num_cols * (num_cols - 1))
		num_rows = num_cols - 1;
	else
		num_rows = num_cols;

	// take the maximum between number of cols and rows
	float max_cols_rows = (num_cols > num_rows) ? num_cols : num_rows;
	// define size of each piece
	float size = 1/ max_cols_rows;

	// get number of column and row for the current index
	int i_col = index % num_cols;
	int i_row = floor(index/num_cols);

	// return the vector already prepared in the format required by glViewport function
	return Vector4(width*i_col*size, height*i_row*size, width*size, height*size);
}

// same function as before but returning the final vector in the format required by the shader
Vector4 GTR::Renderer::assignMapPiece_shader(int width, int height, int index, int num_elements) {
	int num_cols = (int)ceil(sqrt(num_elements));
	int num_rows;
	if (num_elements <= num_cols * (num_cols - 1))
		num_rows = num_cols - 1;
	else
		num_rows = num_cols;

	float max_cols_rows = (num_cols > num_rows) ? num_cols : num_rows;
	float size = 1 / max_cols_rows;

	int i_col = index % num_cols;
	int i_row = floor(index / num_cols);

	return Vector4(size, size, i_col * size, i_row * size);
}

// generate the shadowmap given a light
void GTR::Renderer::generateShadowmap(LightEntity* light, int index)
{
	Camera* view_camera;

	// only spot and directional lights cast shadows
	if (light->light_type != LightEntity::eTypeOfLight::SPOT && light->light_type != LightEntity::eTypeOfLight::DIRECTIONAL)
		return;
	if (!light->cast_shadows) {
		// if a created light is no longer casting shadows delete it
		if (fbo) {
			delete fbo;
			fbo = NULL;
			shadowmap = NULL;
		}
		return;
	}

	if (!fbo) {
		shadowmap = new Texture();
		fbo = new FBO();
		// We only need to store the depth buffer
		width_shadowmap *= num_lights_shadows;
		height_shadowmap *= num_lights_shadows;
		fbo->setDepthOnly(width_shadowmap, height_shadowmap);
		// take the texture from the fbo and store it in another variable
		shadowmap = fbo->depth_texture;
	}

	// Create a new camera from light
	if (!light->light_camera)
		light->light_camera = new Camera();

	// Guardamos la camara anterior para no perderla
	view_camera = Camera::current;
	// activate fbo to start painting in it and not in the screen
	fbo->bind();

	Camera* light_camera = light->light_camera;

	if (light->light_type == LightEntity::eTypeOfLight::SPOT) {
		// set the perspective matrix for the light
		light_camera->setPerspective(light->cone_angle * 2, 1.0, 0.1, light->max_distance);
		// locate and rotate the camera according to the light position, forward direction and up vector
		light_camera->lookAt(light->model.getTranslation(), light->model * Vector3(0, 0, -1), light->model.rotateVector(Vector3(0, 1, 0)));
	}
	else if (light->light_type == LightEntity::eTypeOfLight::DIRECTIONAL) {
		// tried to make the light follow the camera to reach all parts of the scene, but this also makes the light rotate, so it is not working
		//vec3 light_cam_pos = vec3(view_camera->eye.x, view_camera->eye.y, view_camera->eye.z);
		Application* app = Application::instance;
		float halfarea = light->area_size / 2;
		float aspect = Application::instance->window_width / (float)Application::instance->window_height;
		// set orthographic matrix for the light since all rays are parallel
		light_camera->setOrthographic(-halfarea, halfarea, -halfarea * aspect, halfarea * aspect, 0.1, light->max_distance);
		// locate and rotate the camera
		// Now, define center using the target vector since it corresponds to a point where the light is pointing
		light_camera->lookAt(light->model.getTranslation(), light->target, light->model.rotateVector(Vector3(0, 1, 0)));
	}

	light_camera->enable();

	// clear depth buffer to avoid ghosting artifacts only at first light
	if (index == 0)
		glClear(GL_DEPTH_BUFFER_BIT);

	// paint all rendercalls
	for (int i = 0; i < render_calls.size(); i++) {
		RenderCall& rc = render_calls[i];
		// transparent materials do not cast shadows
		if (rc.material->alpha_mode == eAlphaMode::BLEND)
			continue;
		// render if the node is inside the frustum of the new camera
		if (light_camera->testBoxInFrustum(rc.world_bounding.center, rc.world_bounding.halfsize)) {
			// assign position in the shadowmap
			Vector4 piece = assignMapPiece(width_shadowmap, height_shadowmap, index, num_lights_shadows);
			glViewport(piece.x, piece.y, piece.z, piece.w);
			// render shadowmap
			renderFlatMesh(rc.model, rc.mesh, rc.material, light_camera);
		}
	}

	glViewport(0, 0, Application::instance->window_width, Application::instance->window_height);
	fbo->unbind();

	// go back to default system
	view_camera->enable();
	glEnable(GL_DEPTH_TEST);
}

// to show the shadowmap for debugging purposes
void Renderer::showShadowmap(LightEntity* light) {
	if (!shadowmap)
		return;

	Shader* depth_shader = Shader::getDefaultShader("depth");
	depth_shader->enable();
	// set uniforms to delinearize shadowmap texture
	depth_shader->setUniform("u_camera_nearfar", Vector2(light->light_camera->near_plane, light->light_camera->far_plane));
	glViewport(0, 0, 256, 256);

	shadowmap->toViewport(depth_shader);

	// return to default
	glViewport(0, 0, Application::instance->window_width, Application::instance->window_height);
	glEnable(GL_DEPTH_TEST);
}

// --- Render functions ---

// render all entities of the scene
void Renderer::renderScene(GTR::Scene* scene, Camera* camera)
{
	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkGLErrors();

	//render entities
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible)
			continue;

		//is a prefab!
		if (ent->entity_type == PREFAB)
		{
			PrefabEntity* pent = (GTR::PrefabEntity*)ent;
			if (pent->prefab)
				renderPrefab(ent->model, pent->prefab, camera);
		}
	}
}

// To render the scene according to the rendercalls vector
void Renderer::renderScene_RenderCalls(GTR::Scene* scene, Camera* camera) {
	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkGLErrors();

	// Create the lights vector
	lights.clear();
	num_lights_shadows = 0;
	for (int i = 0; i < scene->entities.size(); i++) {
		BaseEntity* ent = scene->entities[i];
		if (ent->entity_type == GTR::eEntityType::LIGHT) {
			LightEntity* light = (LightEntity*)ent;
			if (light->visible) {
				lights.push_back(light);
				if (light->cast_shadows)
					num_lights_shadows += 1;
			}
		}
	}

	// Generate shadowmaps
	// to know if it is the first light casting shadows and therefore clear the depthbuffer of the fbo
	int lshadow_i = 0;
	for (int i = 0; i < lights.size(); i++) {
		LightEntity* light = lights[i];
		// generate only if light is inside the frustum
		if (light->cast_shadows) {
			if (camera->testSphereInFrustum(light->model.getTranslation(), light->max_distance)) {
				generateShadowmap(light, lshadow_i);
				lshadow_i += 1;
			}
		}
	}

	// Create the vector of nodes
	createRenderCalls(scene, camera);

	// Sort the objects by distance to the camera
	sortRenderCalls();

	// check render pipeline
	if (pipeline == FORWARD)
		renderForward(camera);
	else
		renderDeferred(camera);

	if (show_shadowmap)
		showShadowmap(lights[debug_shadowmap]);
}

//renders all the prefab
void Renderer::renderPrefab(const Matrix44& model, GTR::Prefab* prefab, Camera* camera)
{
	assert(prefab && "PREFAB IS NULL");
	//assign the model to the root node
	renderNode(model, &prefab->root, camera);
}

//renders a node of the prefab and its children
void Renderer::renderNode(const Matrix44& prefab_model, GTR::Node* node, Camera* camera)
{
	if (!node->visible)
		return;

	//compute global matrix
	Matrix44 node_model = node->getGlobalMatrix(true) * prefab_model;

	//does this node have a mesh? then we must render it
	if (node->mesh && node->material)
	{
		//compute the bounding box of the object in world space (by using the mesh bounding box transformed to world space)
		BoundingBox world_bounding = transformBoundingBox(node_model, node->mesh->box);

		//if bounding box is inside the camera frustum then the object is probably visible
		if (camera->testBoxInFrustum(world_bounding.center, world_bounding.halfsize))
		{
			//render node mesh
			renderMeshWithMaterial(node_model, node->mesh, node->material, camera);
			//node->mesh->renderBounding(node_model, true);
		}
	}

	//iterate recursively with children
	for (int i = 0; i < node->children.size(); ++i)
		renderNode(prefab_model, node->children[i], camera);
}

//renders a mesh given its transform and material
void Renderer::renderMeshWithMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	Shader* shader = NULL;

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);

	assert(glGetError() == GL_NO_ERROR);

	//chose a shader
	Scene* scene = Scene::instance;
	if (typeOfRender == eRenderPipeline::SINGLEPASS)
		shader = Shader::Get("single_pass");
	else if (typeOfRender == eRenderPipeline::MULTIPASS)
	{
		if (pbr)
			shader = Shader::Get("pbr_multi");
		else
			shader = Shader::Get("multi_pass");
	}

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_model", model);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	// pass textures to the shader
	setTextures(material, shader);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == GTR::eAlphaMode::MASK ? material->alpha_cutoff : 0);

	// pass light parameters
	if (typeOfRender == eRenderPipeline::SINGLEPASS) {
		setSinglepass_parameters(material, shader, mesh);
		//disable shader
		shader->disable();
	}
	else if (typeOfRender == eRenderPipeline::MULTIPASS) {
		setMultipassParameters(material, shader, mesh);
		shader->disable();

	}

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glDepthFunc(GL_LESS);
}

// to render only materials with alpha at the end of deferred rendering
void GTR::Renderer::renderTransparentMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	Shader* shader = NULL;

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);

	assert(glGetError() == GL_NO_ERROR);

	//chose a shader
	if (pbr)
		shader = Shader::Get("pbr_multi");
	else
		shader = Shader::Get("multi_pass");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_model", model);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	// pass textures to the shader
	setTextures(material, shader);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == GTR::eAlphaMode::MASK ? material->alpha_cutoff : 0);

	Scene* scene = Scene::instance;

	Texture* emissive_texture = NULL;
	emissive_texture = material->emissive_texture.texture;

	if (emissive_texture)
		shader->setUniform("u_emissive_texture", emissive_texture, 1);
	// black texture will not add additional light
	else
		shader->setUniform("u_emissive_texture", Texture::getBlackTexture(), 1);

	// paint if value is less or equal to the one in the depth buffer
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	Vector3 ambient_light = scene->ambient_light;
	shader->setUniform("u_ambient_light", ambient_light);
	shader->setUniform("u_light_color", Vector3(0.0, 0.0, 0.0));
	shader->setUniform("u_light_is_first", true);

	mesh->render(GL_TRIANGLES);

	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glDepthFunc(GL_LESS);
}

// to render GBuffers textures
void GTR::Renderer::renderMeshWithMaterialToGBuffers(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (material->alpha_mode == eAlphaMode::BLEND)
		return;

	Shader* shader = NULL;

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	//chose a shader
	Scene* scene = Scene::instance;

	shader = Shader::Get("gbuffers");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_model", model);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);

	Texture* emissive_texture = NULL;
	emissive_texture = material->emissive_texture.texture;

	if (emissive_texture)
		shader->setUniform("u_emissive_texture", emissive_texture, 1);
	// black texture will not add additional light
	else
		shader->setUniform("u_emissive_texture", Texture::getBlackTexture(), 1);

	// pass textures to the shader
	setTextures(material, shader);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == GTR::eAlphaMode::MASK ? material->alpha_cutoff : 0);

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glDepthFunc(GL_LESS);
}

// to save fbo with depth buffer
void Renderer::renderFlatMesh(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera) {
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	Shader* shader = NULL;

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	//chose a shader
	Scene* scene = Scene::instance;
	shader = Shader::Get("flat");


	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_model", model);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == GTR::eAlphaMode::MASK ? material->alpha_cutoff : 0);

	// don't need blending
	glDepthFunc(GL_LESS);
	glDisable(GL_BLEND);

	mesh->render(GL_TRIANGLES);
	//disable shader
	shader->disable();
}

// forward pipeline
void GTR::Renderer::renderForward(Camera* camera)
{
	Scene* scene = Scene::instance;

	// create screen fbo to paint the whole scene and then apply color correction
	if (!screen_fbo) {
		screen_texture = new Texture();
		screen_fbo = new FBO();
		// only need one texture, four channels and depth texture
		screen_fbo->create(Application::instance->window_width, Application::instance->window_height,
			1,
			GL_RGBA,
			GL_UNSIGNED_BYTE,
			true);
	}

	if (typeOfRender == eRenderPipeline::MULTIPASS) {
		// take the texture from the fbo and store it in another variable
		screen_texture = screen_fbo->color_textures[0];

		screen_fbo->bind();
		glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	//render rendercalls
	for (int i = 0; i < render_calls.size(); ++i) {
		// Instead of rendering the entities vector, render the render_calls vector
		RenderCall rc = render_calls[i];

		// if rendercall has mesh and material, render it
		if (rc.mesh && rc.material) {
			// test if node inside the frustum of the camera
			if (camera->testBoxInFrustum(rc.world_bounding.center, rc.world_bounding.halfsize))
				renderMeshWithMaterial(rc.model, rc.mesh, rc.material, camera);
		}
	}

	// gamma correction needs to be applyied after rendering all lights (only in multipass)
	if (typeOfRender == eRenderPipeline::MULTIPASS) {
		screen_fbo->unbind();

		Shader* col_corr = Shader::Get("col_corr");
		col_corr->enable();
		col_corr->setUniform("u_gamma", gamma);
		col_corr->setUniform("u_tonemapping", tonemapping);
		col_corr->setUniform("u_screen_texture", screen_texture, 1);
		col_corr->setUniform("u_lumwhite2", lumwhite2);
		col_corr->setUniform("u_average_lum", averagelum);
		col_corr->setUniform("u_scale", scale);
		col_corr->setUniform("u_depth_texture", screen_fbo->depth_texture, 2);

		Mesh* quad = Mesh::getQuad();
		quad->render(GL_TRIANGLES);
		col_corr->disable();

		glEnable(GL_DEPTH_TEST);
	}
}

// Deferred pipeline
void GTR::Renderer::renderDeferred(Camera* camera)
{
	Scene* scene = Scene::instance;

	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkGLErrors();

	int width = Application::instance->window_width;
	int height = Application::instance->window_height;

	// create FBOs
	if (!gbuffers_fbo) {
		gbuffers_fbo = new FBO();
		illumination_fbo = new FBO();
		ssao_fbo = new FBO();

		//create 3 textures of 4 components. We need more bytes for channel to avoid clipping and losing information
		gbuffers_fbo->create(width, height,
			3,
			GL_RGBA,
			GL_FLOAT,
			true);

		// fbo to compute the illumination for each pixel
		illumination_fbo->create(width, height,
			1,
			GL_RGBA,
			GL_FLOAT,
			true);

		// ssao fbo. RGB to avoid problems
		ssao_fbo->create(width, height,
			1,
			GL_RGB,
			GL_UNSIGNED_BYTE,
			false);
	}

	// generate GBuffers
	gbuffers_fbo->bind();
	generateGBuffers(scene, camera);
	gbuffers_fbo->unbind();

	Matrix44 inv_vp = camera->viewprojection_matrix;
	inv_vp.inverse();

	// generate SSAO
	ssao_fbo->bind();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	generateSSAO(camera, inv_vp, width, height);
	ssao_fbo->unbind();

	// apply illumination
	illumination_fbo->bind();

	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT);

	// don't need to test depth since the objects are already painted
	glDisable(GL_DEPTH_TEST);
	applyIllumination_deferred(scene, camera, inv_vp, width, height);

	// render objects with alpha
	//now we copy the gbuffers depth buffer to the binded one since now we need to test depth
	gbuffers_fbo->depth_texture->copyTo(NULL);
	glEnable(GL_DEPTH_TEST);

	for (int i = 0; i < render_calls.size(); ++i) {
		RenderCall rc = render_calls[i];

		if (rc.mesh && rc.material) {
			if (rc.material->alpha_mode == eAlphaMode::BLEND)
			// test if node inside the frustum of the camera
				if (camera->testBoxInFrustum(rc.world_bounding.center, rc.world_bounding.halfsize))
					renderTransparentMaterial(rc.model, rc.mesh, rc.material, camera);
		}
	}
	glDisable(GL_DEPTH_TEST);
	illumination_fbo->unbind();
	
	glDisable(GL_BLEND);
	applyColorCorrection();
	
	if (show_buffers)
		showGBuffers(Application::instance->window_width, Application::instance->window_height, camera);

	if (show_ssao)
	{
		ssao_fbo->color_textures[0]->toViewport();
	}

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
}

// -- Deferred functions --

// to generate the gbuffers
void GTR::Renderer::generateGBuffers(Scene* scene, Camera* camera)
{
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render in gbuffer
	for (int i = 0; i < render_calls.size(); ++i) {
		// Instead of rendering the entities vector, render the render_calls vector
		RenderCall rc = render_calls[i];

		// if rendercall has mesh and material, render it
		if (rc.mesh && rc.material) {
			// test if node inside the frustum of the camera
			if (camera->testBoxInFrustum(rc.world_bounding.center, rc.world_bounding.halfsize))
				renderMeshWithMaterialToGBuffers(rc.model, rc.mesh, rc.material, camera);
		}
	}
}

// to generate SSAO texture
void GTR::Renderer::generateSSAO(Camera* camera, Matrix44 inv_vp, int width, int height)
{
	Shader* shader = NULL;
	// ssao needs points from a whole sphere
	if (SSAOType == SSAO) {
		shader = Shader::Get("ssao");
		shader->enable();
		shader->setUniform3Array("u_points", (float*)&random_points_sph[0], random_points_sph.size());
	}
	// ssao+ takes points from a hemisphere
	else {
		shader = Shader::Get("ssao_plus");
		shader->enable();
		shader->setUniform3Array("u_points", (float*)&random_points_hemi[0], random_points_hemi.size());

	}
	shader->setUniform("u_gb1_texture", gbuffers_fbo->color_textures[1], 0);
	shader->setUniform("u_depth_texture", gbuffers_fbo->depth_texture, 2);
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_inverse_viewprojection", inv_vp);
	shader->setUniform("u_iRes", Vector2(1.0 / (float)width, 1.0 / (float)height));

	quad->render(GL_TRIANGLES);
	shader->disable();
}

// to apply illumination using deferred
void GTR::Renderer::applyIllumination_deferred(Scene* scene, Camera* camera, Matrix44 inv_vp, int width, int height)
{
	Shader* shader = NULL;
	if (pbr)
		shader = Shader::Get("deferred_pbr");
	else
		shader = Shader::Get("deferred");
	shader->enable();
	shader->setUniform("u_ambient_light", scene->ambient_light);
	shader->setUniform("u_gb0_texture", gbuffers_fbo->color_textures[0], 0);
	shader->setUniform("u_gb1_texture", gbuffers_fbo->color_textures[1], 1);
	shader->setUniform("u_depth_texture", gbuffers_fbo->depth_texture, 4);
	if (ssao)
		shader->setUniform("u_ssao_texture", ssao_fbo->color_textures[0], 5);
	else
		shader->setUniform("u_ssao_texture", Texture::getWhiteTexture(), 5);

	shader->setUniform("u_camera_position", camera->eye);
	//pass the inverse projection of the camera to reconstruct world pos.
	shader->setUniform("u_inverse_viewprojection", inv_vp);
	//pass the inverse window resolution, this may be useful
	shader->setUniform("u_iRes", Vector2(1.0 / (float)width, 1.0 / (float)height));

	shader->setUniform("u_gamma", gamma);

	// temp_ambient will change after first iteration to avoid adding it more than once
	Vector3 temp_ambient = scene->ambient_light;

	glEnable(GL_CULL_FACE);

	bool has_directional = true;
	// if moonlight is not visible, we will not be rendering all the screen, only the spheres where lights are
	// so we need to render first the quad with only ambient and then add the rest of lights
	if (lights.size()) {
		// check if directional light is visible
		int lights_size = lights.size();
		has_directional = lights[lights_size - 1]->light_type == LightEntity::eTypeOfLight::DIRECTIONAL && lights[lights_size - 1]->visible;
	}
	// render a quad with only ambient light if there are no lights or no directional light
	if (!lights.size() || !has_directional) {
		shader->disable();
		// Creo otro shader por qu?? ahora pintaremos un quad en lugar de esferas, y entonces necesitaremos el quad.vs y las uv habituales de la textura
		Shader* shader_ambient = Shader::Get("deferred_ambient");
		shader_ambient->enable();
		shader_ambient->setUniform("u_ambient_light", scene->ambient_light);
		shader_ambient->setUniform("u_gamma", gamma);
		shader_ambient->setUniform("u_gb0_texture", gbuffers_fbo->color_textures[0], 0);
		shader_ambient->setUniform("u_gb1_texture", gbuffers_fbo->color_textures[1], 1);
		shader_ambient->setUniform("u_gb2_texture", gbuffers_fbo->color_textures[2], 2);
		if (ssao)
			shader->setUniform("u_ssao_texture", ssao_fbo->color_textures[0], 5);
		else
			shader->setUniform("u_ssao_texture", Texture::getWhiteTexture(), 5);
		quad->render(GL_TRIANGLES);
		shader_ambient->disable();
		// set ambient light to zero to avoid adding it again
		temp_ambient = Vector3(0.0, 0.0, 0.0);
	}
	shader->enable();
	int i_shadow = 0;

	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_GREATER);
	////block writing to the ZBuffer so we do not modify it with our geometry
	//glDepthMask(false);

	// render the rest of lights
	for (int i = 0; i < lights.size(); i++) {
		// emissive texture -- add it only once
		if (i == 0)
			shader->setUniform("u_gb2_texture", gbuffers_fbo->color_textures[2], 2);
		else
		{
			shader->setUniform("u_gb2_texture", Texture::getBlackTexture(), 2);
		}

		LightEntity* light = lights[i];

		Vector3 dist_cam_light = light->model.getTranslation() - camera->center;
		float module_dist = dist_cam_light.length();

		// check which side of the sphere we need to render
		if (module_dist < light->max_distance)
			glFrontFace(GL_CW);
		else
			glFrontFace(GL_CCW);

		if (i == 0) {
			glDisable(GL_BLEND);
		}
		else {
			glBlendFunc(GL_ONE, GL_ONE);
			glEnable(GL_BLEND);
		}
		sphere->radius = light->max_distance;

		// set position and scale of the sphere according to light
		Matrix44 m;
		Vector3 light_pos = light->model.getTranslation();
		m.setTranslation(light_pos.x, light_pos.y, light_pos.z);
		m.scale(light->max_distance, light->max_distance, light->max_distance);
		shader->setUniform("u_model", m);
		shader->setUniform("u_viewprojection", camera->viewprojection_matrix);

		uploadLightToShader(light, shader, temp_ambient, i_shadow);
		shader->setUniform("u_light_is_last", i == lights.size() - 1 ? 1 : 0);

		// increment the number of lights casting shadows
		if (light->cast_shadows) {
			i_shadow += 1;
		}

		sphere->render(GL_TRIANGLES);

		// change ambient after first iteration
		temp_ambient = Vector3(0.0, 0.0, 0.0);
		glFrontFace(GL_CCW);

	}
	shader->disable();
	//glDisable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);
	////block writing to the ZBuffer so we do not modify it with our geometry
	//glDepthMask(true);
	//glFrontFace(GL_CCW);
}

// to apply color correction given the scene texture
void GTR::Renderer::applyColorCorrection()
{
	// apply color correction
	Shader* col_corr = Shader::Get("col_corr");
	col_corr->enable();
	col_corr->setUniform("u_gamma", gamma);
	col_corr->setUniform("u_screen_texture", illumination_fbo->color_textures[0], 0);
	col_corr->setUniform("u_tonemapping", tonemapping);
	col_corr->setUniform("u_lumwhite2", lumwhite2);
	col_corr->setUniform("u_average_lum", averagelum);
	col_corr->setUniform("u_scale", scale);
	col_corr->setUniform("u_depth_texture", illumination_fbo->depth_texture, 1);

	Mesh* quad_final = Mesh::getQuad();
	quad_final->render(GL_TRIANGLES);
	col_corr->disable();
}

// -- Upload to shader functions --

void Renderer::uploadLightToShader(GTR::LightEntity* light, Shader* shader, Vector3 ambient_light, int shadow_i) {	
	shader->setUniform("u_light_type", light->light_type);
	shader->setUniform("u_ambient_light", ambient_light);
	shader->setUniform("u_light_position", light->model.getTranslation());
	shader->setUniform("u_light_color", light->color * light->intensity);
	shader->setUniform("u_light_max_distance", light->max_distance);

	// Use the cosine to compare it directly to NdotL
	shader->setUniform("u_light_cone_cos", (float)cos(light->cone_angle * DEG2RAD));
	shader->setUniform("u_light_cone_exp", light->cone_exp);

	if (light->light_type == LightEntity::eTypeOfLight::SPOT)
		shader->setUniform("u_light_direction", light->model.rotateVector(Vector3(0.0, 0.0, -1.0)));
	else if (light->light_type == LightEntity::eTypeOfLight::DIRECTIONAL)
		shader->setUniform("u_light_direction", (light->model.getTranslation() - light->target));

	if (shadowmap && light->cast_shadows) {
		shader->setUniform("u_light_cast_shadows", 1);
		shader->setUniform("u_light_shadowmap", shadowmap, 8);
		shader->setUniform("u_light_shadowmap_vpm", light->light_camera->viewprojection_matrix);
		shader->setUniform("u_light_shadow_bias", light->shadow_bias);
		shader->setUniform("u_shadow_i", shadow_i);
		vec4 map_piece = assignMapPiece_shader(2048 * 3, 2048 * 3, shadow_i, num_lights_shadows);
		shader->setUniform("u_map_piece", map_piece);
	}
	else
		shader->setUniform("u_light_cast_shadows", 0);
}

// to upload textures to the shader
void Renderer::setTextures(GTR::Material* material, Shader* shader) {
	Texture* texture = NULL;
	Texture* normal_texture = NULL;
	Texture* occlusion_texture = NULL;
	Texture* met_rough_texture = NULL;

	// save textures
	texture = material->color_texture.texture;
	normal_texture = material->normal_texture.texture;
	met_rough_texture = material->metallic_roughness_texture.texture;
	occlusion_texture = material->occlusion_texture.texture;

	if (texture == NULL)
		texture = Texture::getWhiteTexture(); //a 1x1 white texture

	// pass textures
	if (texture)
		shader->setUniform("u_texture", texture, 0);

	if (occlusion_texture)
		shader->setUniform("u_occlusion_texture", occlusion_texture, 2);
	// white texture will take into account all light, namely no occlusion
	else
		shader->setUniform("u_occlusion_texture", Texture::getWhiteTexture(), 2);

	if (met_rough_texture)
		shader->setUniform("u_met_rough_texture", met_rough_texture, 3);
	else
		shader->setUniform("u_met_rough_texture", Texture::getWhiteTexture(), 3);

	if (normal_texture) {
		shader->setUniform("u_normal_texture", normal_texture, 4);
		shader->setUniform("u_normal_text_bool", 1);
	}
	// if the material do not have normal texture as the floor, we set the boolean to false to avoid artifacts
	else
		shader->setUniform("u_normal_text_bool", 0);

	shader->setUniform("u_texture2show", debug_texture);
}

void Renderer::setSinglepass_parameters(GTR::Material* material, Shader* shader, Mesh* mesh) {
	Texture* emissive_texture = NULL;
	emissive_texture = material->emissive_texture.texture;

	//select the blending
	if (material->alpha_mode == GTR::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	if (emissive_texture)
		shader->setUniform("u_emissive_texture", emissive_texture, 1);
	// black texture will not add additional light
	else
		shader->setUniform("u_emissive_texture", Texture::getBlackTexture(), 1);

	// Define some variables to store lights information
	std::vector<int> lights_type;
	std::vector<Vector3> lights_position;
	std::vector<Vector3> lights_color;
	std::vector<float> lights_max_distance;
	// spot and directional parameters
	std::vector<float> lights_cone_cos;
	std::vector<float> lights_cone_exp;
	std::vector<Vector3> lights_direction;
	// shadows parameters
	std::vector<int> lights_cast_shadows;
	// didn't make it to pass multiple textures to the shader, so shadows are not working in singlepass
	//std::vector<Texture*> lights_shadowmap;
	std::vector<Matrix44> lights_shadowmap_vpm;
	std::vector<float> lights_shadow_bias;

	int lshadow_count = 0;
	// Iterate and store the information
	for (int i = 0; i < max_lights; i++) {
		LightEntity* light;
		// if the light is not created, pass an empty light to the shader to avoid taking trash information
		int lights_size = 0;
		if (lights.size()) {
			lights_size = lights.size();
		}
		if (i < lights_size && lights[i]->visible)
			light = lights[i];
		else
			light = new LightEntity();

		// add the information to the vectors
		lights_type.push_back(light->light_type);
		lights_position.push_back(light->model.getTranslation());
		lights_color.push_back(light->color * light->intensity);
		lights_max_distance.push_back(light->max_distance);
		lights_cone_cos.push_back((float)cos(light->cone_angle * DEG2RAD));
		lights_cone_exp.push_back(light->cone_exp);

		if (light->light_type == LightEntity::eTypeOfLight::SPOT)
			// take the forward vector of the light
			lights_direction.push_back(light->model.rotateVector(Vector3(0.0, 0.0, -1.0)));
		else if (light->light_type == LightEntity::eTypeOfLight::DIRECTIONAL)
			// get direction of the light taking the target point
			lights_direction.push_back((light->model.getTranslation() - light->target));
		else {
			lights_direction.push_back(vec3(0.0, 0.0, 0.0));
		}

		lights_shadow_bias.push_back(light->shadow_bias);
		if (light->cast_shadows) {
			lights_cast_shadows.push_back(1);
			//lights_shadowmap.push_back(light->shadowmap);
			lights_shadowmap_vpm.push_back(light->light_camera->viewprojection_matrix);
		}
		else {
			lights_cast_shadows.push_back(0);
		}
	}

	// Pass to the shader
	Scene* scene = Scene::instance;

	shader->setUniform("u_gamma", gamma);

	if (!lights.size()) {
		int t = 0;
	}
	shader->setUniform("u_lights_type", lights_type);
	shader->setUniform("u_ambient_light", scene->ambient_light);
	shader->setUniform("u_lights_position", lights_position);
	shader->setUniform("u_lights_color", lights_color);
	shader->setUniform("u_lights_max_distance", lights_max_distance);

	// Use the cosine to compare it directly to NdotL
	shader->setUniform("u_lights_cone_cos", lights_cone_cos);
	shader->setUniform("u_lights_cone_exp", lights_cone_exp);
	shader->setUniform("u_lights_direction", lights_direction);

	shader->setUniform("u_num_lights_shadow", num_lights_shadows);
	shader->setUniform("u_lights_cast_shadows", lights_cast_shadows);
	shader->setUniform("u_lights_shadowmap", shadowmap, 8);
	if (lights_shadowmap_vpm.size())
		shader->setUniform("u_lights_shadowmap_vpm", lights_shadowmap_vpm);
	shader->setUniform("u_lights_shadow_bias", lights_shadow_bias);

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	// clear all vectors
	lights_type.clear();
	lights_position.clear();
	lights_color.clear();
	lights_max_distance.clear();
	lights_cone_cos.clear();
	lights_cone_exp.clear();
	lights_direction.clear();
	lights_cast_shadows.clear();
	lights_shadowmap_vpm.clear();
	lights_shadow_bias.clear();
}

void Renderer::setMultipassParameters(GTR::Material* material, Shader* shader, Mesh* mesh) {
	Scene* scene = Scene::instance;


	Texture* emissive_texture = NULL;
	emissive_texture = material->emissive_texture.texture;

	if (emissive_texture)
		shader->setUniform("u_emissive_texture", emissive_texture, 1);
	// black texture will not add additional light
	else
		shader->setUniform("u_emissive_texture", Texture::getBlackTexture(), 1);

	// paint if value is less or equal to the one in the depth buffer
	glDepthFunc(GL_LEQUAL);

	// do a linear interpolation adding the pixel painted with the current one
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	Vector3 temp_ambient = scene->ambient_light;

	// to know if we are in first iteration of a visible light, since the first light can be disabled and therefore blending will not be activated for transparent materials
	bool is_first = true;

	if (!lights.size()) {
		shader->setUniform("u_ambient_light", temp_ambient);
		shader->setUniform("u_light_color", Vector3(0.0, 0.0, 0.0));

		mesh->render(GL_TRIANGLES);
	}

	int lshadow_count = 0;
	for (int i = 0; i < lights.size(); ++i) {
		LightEntity* light = lights[i];

		if (!lights[i]->visible)
			continue;

		if (is_first) {
			// select the blending
			if (material->alpha_mode == GTR::eAlphaMode::BLEND)
			{
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			}
			else
				glDisable(GL_BLEND);
		}

		// if a texture is selected in imgui, render just one time
		else if (debug_texture != eTextureType::COMPLETE)
			continue;

		// we already passed first light

		shader->setUniform("u_gamma", gamma);
		shader->setUniform("u_light_is_first", is_first);

		// Pass to the shader
		shader->setUniform("u_light_type", light->light_type);
		shader->setUniform("u_ambient_light", temp_ambient);
		shader->setUniform("u_light_position", light->model.getTranslation());
		shader->setUniform("u_light_color", light->color * light->intensity);
		shader->setUniform("u_light_max_distance", light->max_distance);

		// Use the cosine to compare it directly to NdotL
		shader->setUniform("u_light_cone_cos", (float)cos(light->cone_angle * DEG2RAD));
		shader->setUniform("u_light_cone_exp", light->cone_exp);

		if (light->light_type == LightEntity::eTypeOfLight::SPOT)
			shader->setUniform("u_light_direction", light->model.rotateVector(Vector3(0.0, 0.0, -1.0)));
		else if (light->light_type == LightEntity::eTypeOfLight::DIRECTIONAL)
			shader->setUniform("u_light_direction", (light->model.getTranslation() - light->target));

		if (light->cast_shadows) {
			shader->setUniform("u_light_cast_shadows", 1);
			shader->setUniform("u_light_shadowmap", shadowmap, 8);
			shader->setUniform("u_light_shadowmap_vpm", light->light_camera->viewprojection_matrix);
			shader->setUniform("u_light_shadow_bias", light->shadow_bias);
			shader->setUniform("u_shadow_i", lshadow_count);
			vec4 map_piece = assignMapPiece_shader(width_shadowmap, height_shadowmap, lshadow_count, num_lights_shadows);
			shader->setUniform("u_map_piece", map_piece);
			lshadow_count += 1;
		}
		else
			shader->setUniform("u_light_cast_shadows", 0);

		is_first = false;

		//do the draw call that renders the mesh into the screen
		mesh->render(GL_TRIANGLES);

		// Activate blending again for the rest of lights to do the interpolation
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE);

		// Reset ambient light to add it only once
		temp_ambient = vec3(0.0, 0.0, 0.0);
		emissive_texture = Texture::getBlackTexture();
		shader->setUniform("u_emissive_texture", emissive_texture, 1);
	}
}

// -- Debug functions --

void Renderer::showGBuffers(int width, int height, Camera* camera) {
	// to avoid seeing the scene
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// four textures
	float size = 0.5;

	glViewport(0, height * size, width * size, height * size);
	Shader* depth_shader = Shader::getDefaultShader("depth");
	depth_shader->enable();
	// set uniforms to delinearize shadowmap texture
	depth_shader->setUniform("u_camera_nearfar", Vector2(camera->near_plane, camera->far_plane));
	gbuffers_fbo->depth_texture->toViewport(depth_shader);

	glViewport(width * size, height * size, width * size, height * size);
	gbuffers_fbo->color_textures[0]->toViewport();

	glViewport(0, 0, width * size, height * size);
	gbuffers_fbo->color_textures[1]->toViewport();

	glViewport(width * size, 0, width * size, height * size);
	gbuffers_fbo->color_textures[2]->toViewport();

	glViewport(0, 0, width, height);
	glEnable(GL_DEPTH_TEST);
}

// to set all lights to visible
void GTR::Renderer::setLightsVisible()
{
	Scene* scene = Scene::instance;
	for (int i = 0; i < scene->entities.size(); i++) {
		BaseEntity* ent = scene->entities[i];
		if (ent->entity_type == GTR::eEntityType::LIGHT) {
			LightEntity* light = (LightEntity*)ent;
			light->visible = true;
		}
	}
}

// to set all lights to not visible
void GTR::Renderer::setLightsInvisible()
{
	for (int i = 0; i < lights.size(); i++) {
		lights[i]->visible = false;
	}
}

void GTR::Renderer::renderInMenu() {
	if (ImGui::Button("Bona Nit"))
		setLightsInvisible();
	if (ImGui::Button("Bon Dia"))
		setLightsVisible();
	ImGui::Checkbox("Show Shadowmap", &show_shadowmap);
	if (pipeline == ePipeline::DEFERRED) {
		ImGui::Checkbox("Show GBuffers", &show_buffers);
	}
	else {
		ImGui::Combo("Textures", &debug_texture, "COMPLETE\0NORMAL\0OCCLUSION\0EMISSIVE\0METALNESS\0ROUGHNESS");
	}
}

Texture* GTR::CubemapFromHDRE(const char* filename)
{
	HDRE* hdre = HDRE::Get(filename);
	if (!hdre)
		return NULL;

	Texture* texture = new Texture();
	if (hdre->getFacesf(0))
	{
		texture->createCubemap(hdre->width, hdre->height, (Uint8**)hdre->getFacesf(0),
			hdre->header.numChannels == 3 ? GL_RGB : GL_RGBA, GL_FLOAT);
		for (int i = 1; i < hdre->levels; ++i)
			texture->uploadCubemap(texture->format, texture->type, false,
				(Uint8**)hdre->getFacesf(i), GL_RGBA32F, i);
	}
	else
		if (hdre->getFacesh(0))
		{
			texture->createCubemap(hdre->width, hdre->height, (Uint8**)hdre->getFacesh(0),
				hdre->header.numChannels == 3 ? GL_RGB : GL_RGBA, GL_HALF_FLOAT);
			for (int i = 1; i < hdre->levels; ++i)
				texture->uploadCubemap(texture->format, texture->type, false,
					(Uint8**)hdre->getFacesh(i), GL_RGBA16F, i);
		}
	return texture;
}

std::vector<Vector3> GTR::generateSpherePoints(int num, float radius, bool hemi)
{
	std::vector<Vector3> points;
	points.resize(num);
	for (int i = 0; i < num; i += 1)
	{
		Vector3& p = points[i];
		float u = random();
		float v = random();
		float theta = u * 2.0 * PI;
		float phi = acos(2.0 * v - 1.0);
		float r = cbrt(random() * 0.9 + 0.1) * radius;
		float sinTheta = sin(theta);
		float cosTheta = cos(theta);
		float sinPhi = sin(phi);
		float cosPhi = cos(phi);
		p.x = r * sinPhi * cosTheta;
		p.y = r * sinPhi * sinTheta;
		p.z = r * cosPhi;
		if (hemi && p.z < 0)
			p.z *= -1.0;
	}
	return points;
}

void Renderer::renderProbe(Vector3 pos, float size, float* coeffs)
{
	Camera* camera = Camera::current;
	Shader* shader = Shader::Get("probe");
	Mesh* mesh = Mesh::Get("data/meshes/sphere.obj");

	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	Matrix44 model;
	model.setTranslation(pos.x, pos.y, pos.z);
	model.scale(size, size, size);

	shader->enable();
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_model", model);
	shader->setUniform3Array("u_coeffs", coeffs, 9);

	mesh->render(GL_TRIANGLES);
}
