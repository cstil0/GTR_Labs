#pragma once
#include "prefab.h"
#include "shader.h"
#include <string>

//forward declarations
class Camera;

namespace GTR {

	class Prefab;
	class Material;

	class RenderCall {
	public:
		Mesh* mesh;
		Material* material;
		Matrix44 model;
		float distance_to_camera;
		BoundingBox world_bounding;

		RenderCall() {}
		virtual ~RenderCall() {}
	};

	
	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
		enum eTextureType {
			COMPLETE,
			NORMAL,
			OCCLUSION,
			EMISSIVE
		};

		enum ePipeline {
			FORWARD,
			DEFERRED
		};

	public:

		// Save only the visible nodes sorted by distance to the camera
		std::vector<RenderCall> render_calls;
		// Save all lights in the scene
		std::vector<LightEntity*> lights;

		int max_lights;

		FBO* fbo;
		Texture* shadowmap;
		ePipeline pipeline;

		// Imgui debug parameters
		bool show_shadowmap;
		int debug_shadowmap;
		int debug_texture;

		FBO* gbuffers_fbo;
		FBO* illumination_fbo;


		Renderer();

		// -- Rendercalls manager functions--
		void createRenderCalls(GTR::Scene* scene, Camera* camera);
		void addRenderCall_node(GTR::Scene* scene, Camera* camera, Node* node, Matrix44 curr_model, Matrix44 parent_model);
		void sortRenderCalls();
		// operator used to sort rendercalls vector
		static bool compare_distances(const RenderCall rc1, const RenderCall rc2) { return (rc1.distance_to_camera < rc2.distance_to_camera); }

		// -- Shadowmap functions --
		void showShadowmap(LightEntity* light);
		void generateShadowmap(LightEntity* light);

		// -- Render functions --
		//renders several elements of the scene
		void renderScene(GTR::Scene* scene, Camera* camera);
		// to render the scene using rendercalls vector
		void renderScene_RenderCalls(GTR::Scene* scene, Camera* camera);
		//to render a whole prefab (with all its nodes)
		void renderPrefab(const Matrix44& model, GTR::Prefab* prefab, Camera* camera);
		//to render one node from the prefab and its children
		void renderNode(const Matrix44& model, GTR::Node* node, Camera* camera);
		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		void renderMeshWithMaterialToGBuffers(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		
		void setTextures(GTR::Material* material, Shader* shader);
		void setSinglepass_parameters(GTR::Material* material, Shader* shader, Mesh* mesh);
		void setMultipassParameters(GTR::Material* material, Shader* shader, Mesh* mesh);
		// to render flat objects for generating the shadowmaps
		void renderFlatMesh(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		void showGBuffers(int width, int height, Camera* camera);
		void renderForward(Camera* camera);
		void renderDeferred(Camera* camera);
		void renderInMenu();
	};

	Texture* CubemapFromHDRE(const char* filename);

};