#ifndef SCENE_H
#define SCENE_H

#include "framework.h"
#include "camera.h"
#include "material.h"
#include <string>
#include <map>

//forward declaration
class cJSON; 
class FBO;
class Texture;

//our namespace
namespace GTR {
	enum eEntityType {
		NONE = 0,
		PREFAB = 1,
		LIGHT = 2,
		CAMERA = 3,
		REFLECTION_PROBE = 4,
		DECAL = 5,
	};

	class Scene;
	class Prefab;
	//class Light;
	class Node;

	//represents one element of the scene (could be lights, prefabs, cameras, etc)
	class BaseEntity
	{
	public:
		Scene* scene;
		std::string name;
		eEntityType entity_type;
		Matrix44 model;
		bool visible;

		BaseEntity() { entity_type = NONE; visible = true; }
		virtual ~BaseEntity() {}

		virtual void renderInMenu();
		virtual void configure(cJSON* json) {}
	};

	//represents one prefab in the scene
	class PrefabEntity : public GTR::BaseEntity
	{
	public:
		std::string filename;
		Prefab* prefab;

		PrefabEntity();
		virtual void renderInMenu();
		virtual void configure(cJSON* json);
	};

	// represent a light in the scene
	class LightEntity : public GTR::BaseEntity {
	public:
		enum eTypeOfLight {
			NONE,
			POINT,
			SPOT,
			DIRECTIONAL
		};

		std::string filename;
		Vector3 color;
		float intensity;
		int light_type;
		//float angle;
		bool volumetric;

		int shadow_i;
		bool cast_shadows;
		float shadow_bias;

		float max_distance;
		float cone_angle;
		float cone_exp;
		float area_size;
		Vector3 target;

		FBO* fbo;
		Texture* shadowmap;
		Camera* light_camera;

		LightEntity();

		virtual void renderInMenu();
		virtual void configure(cJSON* json);
	};

	class DecalEntity : public BaseEntity {
	public:
		std::string texture;

		DecalEntity();
		virtual void renderInMenu() {};
		virtual void configure(cJSON* json);
	};

	//class ReflectionProbeEntity : public GTR::BaseEntity {
	//public:
	//	Texture* texture;

	//	ReflectionProbeEntity();
	//	virtual void renderInMenu();
	//	virtual void configure(cJSON* json);
	//};

	//contains all entities of the scene
	class Scene
	{
	public:
		enum eSceneType {
			DEFAULT,
			SPONZA
		};

		static Scene* instance;

		Vector3 background_color;
		Vector3 ambient_light;
		float air_density;
		Camera main_camera;
		int scene_type;
		Texture* skybox;

		std::string filename;
		// entities in the scene
		std::vector<BaseEntity*> entities;
		std::map<std::string, BaseEntity*> entities_by_name;

		Scene();

		void clear();
		void addEntity(BaseEntity* entity);

		BaseEntity* getEntityByName(std::string name);

		bool load(const char* filename);
		BaseEntity* createEntity(std::string type);
	};
	Texture* CubemapFromHDRE(const char* filename);

};

#endif