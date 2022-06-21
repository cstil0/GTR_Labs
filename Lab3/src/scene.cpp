#include "scene.h"
#include "utils.h"

#include "prefab.h"
#include "extra/cJSON.h"
  

GTR::Scene* GTR::Scene::instance = NULL;

GTR::Scene::Scene()
{
	instance = this;
	air_density = 1.0;
	// Start with singlepass
}

void GTR::Scene::clear()
{
	for (int i = 0; i < entities.size(); ++i)
	{
		BaseEntity* ent = entities[i];
		delete ent;
	}
	entities.resize(0);
}

void GTR::Scene::addEntity(BaseEntity* entity)
{
	entities.push_back(entity); 
	entity->scene = this;
	if (entity->name.size())
		entities_by_name[entity->name] = entity;
}

GTR::BaseEntity* GTR::Scene::getEntityByName(std::string name)
{
	auto it = entities_by_name.find(name);
	if (it == entities_by_name.end())
		return NULL;
	return it->second;
}

bool GTR::Scene::load(const char* filename)
{
	std::string content;

	this->filename = filename;
	std::cout << " + Reading scene JSON: " << filename << "..." << std::endl;

	if (!readFile(filename, content))
	{
		std::cout << "- ERROR: Scene file not found: " << filename << std::endl;
		return false;
	}

	//parse json string 
	cJSON* json = cJSON_Parse(content.c_str());
	if (!json)
	{
		std::cout << "ERROR: Scene JSON has errors: " << filename << std::endl;
		return false;
	}

	//read global properties
	background_color = readJSONVector3(json, "background_color", background_color);
	ambient_light = readJSONVector3(json, "ambient_light", ambient_light );
	main_camera.eye = readJSONVector3(json, "camera_position", main_camera.eye);
	main_camera.center = readJSONVector3(json, "camera_target", main_camera.center);
	main_camera.fov = readJSONNumber(json, "camera_fov", main_camera.fov);
	std::string skybox_filename = cJSON_GetObjectItem(json, "environment")->valuestring;
	skybox = CubemapFromHDRE((std::string("data/") + skybox_filename).c_str());
	if (!skybox)
		std::cout << "Skybox texture not found: " << filename << std::endl;

	//entities
	cJSON* entities_json = cJSON_GetObjectItemCaseSensitive(json, "entities");
	cJSON* entity_json;
	cJSON_ArrayForEach(entity_json, entities_json)
	{
		std::string type_str = cJSON_GetObjectItem(entity_json, "type")->valuestring;
		BaseEntity* ent = createEntity(type_str);
		if (!ent)
		{
			std::cout << " - ENTITY TYPE UNKNOWN: " << type_str << std::endl;
			//continue;
			ent = new BaseEntity();
		}

		addEntity(ent);

		if (cJSON_GetObjectItem(entity_json, "name"))
		{
			ent->name = cJSON_GetObjectItem(entity_json, "name")->valuestring;
			stdlog(std::string(" + entity: ") + ent->name);
		}

		//read transform
		if (cJSON_GetObjectItem(entity_json, "position"))
		{
			ent->model.setIdentity();
			Vector3 position = readJSONVector3(entity_json, "position", Vector3());
			ent->model.translate(position.x, position.y, position.z);
		}

		// Add rotation for x and y axis
		if (cJSON_GetObjectItem(entity_json, "angle_x"))
		{
			float angle = cJSON_GetObjectItem(entity_json, "angle_x")->valuedouble;
			ent->model.rotate(angle * DEG2RAD, Vector3(1, 0, 0));
		}

		if (cJSON_GetObjectItem(entity_json, "angle_y"))
		{
			float angle = cJSON_GetObjectItem(entity_json, "angle_y")->valuedouble;
			ent->model.rotate(angle * DEG2RAD, Vector3(0, 0, -1));
		}

		// z axis rotation
		if (cJSON_GetObjectItem(entity_json, "angle"))
		{
			float angle = cJSON_GetObjectItem(entity_json, "angle")->valuedouble;
			ent->model.rotate(angle * DEG2RAD, Vector3(0, 1, 0));
		}

		if (cJSON_GetObjectItem(entity_json, "rotation"))
		{
			Vector4 rotation = readJSONVector4(entity_json, "rotation");
			Quaternion q(rotation.x, rotation.y, rotation.z, rotation.w);
			Matrix44 R;
			q.toMatrix(R);
			ent->model = R * ent->model;
		}

		if (cJSON_GetObjectItem(entity_json, "target"))
		{
			Vector3 target = readJSONVector3(entity_json, "target", Vector3());
			Vector3 front = target - ent->model.getTranslation();
			ent->model.setFrontAndOrthonormalize(front);
		}

		if (cJSON_GetObjectItem(entity_json, "scale"))
		{
			Vector3 scale = readJSONVector3(entity_json, "scale", Vector3(1, 1, 1));
			ent->model.scale(scale.x, scale.y, scale.z);
		}

		ent->configure(entity_json);
	}

	//free memory
	cJSON_Delete(json);

	return true;
}

GTR::BaseEntity* GTR::Scene::createEntity(std::string type)
{
	if (type == "PREFAB")
		return new GTR::PrefabEntity();

	if (type == "LIGHT")
		return new GTR::LightEntity();

	if (type == "DECAL")
		return new GTR::DecalEntity();

	return NULL;
}


void GTR::BaseEntity::renderInMenu()
{
#ifndef SKIP_IMGUI
	ImGui::Text("Name: %s", name.c_str()); // Edit 3 floats representing a color
	ImGui::Checkbox("Visible", &visible); // Edit 3 floats representing a color
	//Model edit
	ImGuiMatrix44(model, "Model");
#endif
}

GTR::PrefabEntity::PrefabEntity()
{
	entity_type = PREFAB;
	prefab = NULL;
}

void GTR::PrefabEntity::configure(cJSON* json)
{
	if (cJSON_GetObjectItem(json, "filename"))
	{
		filename = cJSON_GetObjectItem(json, "filename")->valuestring;
		prefab = GTR::Prefab::Get( (std::string("data/") + filename).c_str());
	}
}

void GTR::PrefabEntity::renderInMenu()
{
	BaseEntity::renderInMenu();

#ifndef SKIP_IMGUI
	ImGui::Text("filename: %s", filename.c_str()); // Edit 3 floats representing a color
	if (prefab && ImGui::TreeNode(prefab, "Prefab Info"))
	{
		prefab->root.renderInMenu();
		ImGui::TreePop();
	}
#endif
}

GTR::LightEntity::LightEntity()
{
	// Initialize to zero all variables
	entity_type = LIGHT;
	light_type = LightEntity::eTypeOfLight::NONE;
	color.set(0, 0, 0);
	intensity = 0;
	max_distance = 0;
	cone_angle = 0;
	cone_exp = 0;

	area_size = 0;
	target = vec3(0.0, 0.0, 0.0);

	shadow_i = -1;
	cast_shadows = false;
	shadow_bias = 0;

	light_camera = NULL;
	fbo = NULL;
	shadowmap = NULL;
	light_camera = NULL;

	volumetric = false;
}

void GTR::LightEntity::renderInMenu() {
	//First render the base menu
	BaseEntity::renderInMenu();
#ifndef SKIP_IMGUI
	ImGui::Text("filename: %s", filename.c_str()); // Edit 3 floats representing a color

	std::string type_str;
	switch (light_type) {
	case eTypeOfLight::POINT: type_str = "POINT"; break;
	case eTypeOfLight::SPOT: type_str = "SPOT"; break;
	case eTypeOfLight::DIRECTIONAL: type_str = "DIRECTIONAL"; break;
	}

	ImGui::Text("LightType: %s", type_str.c_str());
	ImGui::ColorEdit3("Color", color.v);
	ImGui::SliderFloat("Intensity", &intensity, 0.0, 10);
	ImGui::SliderFloat("Maximum Distance", &max_distance, 0.0, 1000);

	// Show the parameters depending on the type of light
	if (light_type == LightEntity::eTypeOfLight::SPOT) {
		ImGui::SliderFloat("Cone Angle", &cone_angle, 0.0, 80);
		ImGui::SliderFloat("Cone Exponential", &cone_exp, 0.0, 100);
		ImGui::Checkbox("Shadows", &cast_shadows);
		ImGui::Checkbox("Volumetric", &volumetric);

	}
	else if (light_type == LightEntity::eTypeOfLight::DIRECTIONAL){
		ImGui::SliderFloat("Area Size", &area_size, 0.0, 2000);
		ImGui::DragFloat3("Target", &target.x, 1, -80, 80);
		ImGui::Checkbox("Shadows", &cast_shadows);
		ImGui::Checkbox("Volumetric", &volumetric);
	}
	if (cast_shadows) {
		ImGui::SliderFloat("Shadow Bias", &shadow_bias, 0.001, 0.5);
	}
#endif

}

void GTR::LightEntity::configure(cJSON* json)
{
	// Read parameters
	color = readJSONVector3(json, "color", color);
	intensity = readJSONNumber(json, "intensity", intensity);
	std::string str = readJSONString(json, "light_type", "");
	if (str == "POINT")
		light_type = eTypeOfLight::POINT;
	else if (str == "SPOT")
		light_type = eTypeOfLight::SPOT;
	else if (str == "DIRECTIONAL")
		light_type = eTypeOfLight::DIRECTIONAL;
	else
		light_type = eTypeOfLight::NONE;

	//angle = readJSONNumber(json, "angle", angle);
	area_size = readJSONNumber(json, "area_size", area_size);
	max_distance = readJSONNumber(json, "max_dist", max_distance);
	cone_angle = readJSONNumber(json, "cone_angle", cone_angle);
	cone_exp = readJSONNumber(json, "cone_exp", cone_angle);
	cast_shadows = readJSONBool(json, "cast_shadows", false);
	shadow_bias = readJSONNumber(json, "shadow_bias", shadow_bias);
	volumetric = readJSONBool(json, "volumetric", volumetric);
}

GTR::DecalEntity::DecalEntity()
{
	entity_type = DECAL;
}

void GTR::DecalEntity::configure(cJSON* json)
{
	if (cJSON_GetObjectItem(json, "texture"))
		texture = cJSON_GetObjectItem(json, "texture")->valuestring;
}
