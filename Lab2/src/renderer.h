#pragma once
#include "prefab.h"
#include "shader.h"
#include "sphericalharmonics.h"
#include <string>

//forward declarations
class Camera;

namespace GTR {

	class Prefab;
	class Material;


	//the struct that holds one probe coeffs
	//already defined in "sphericalharmonics.h"
	// struct SphericalHarmonics {
	//	Vector3 coeffs[9];
	// };

	//struct to store probes
	struct sProbe {
		Vector3 pos; //where is located
		Vector3 local; //its ijk pos in the matrix
		int index; //its index in the linear array
		SphericalHarmonics sh; //coeffs
	};


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
	public:
		// Enums
		enum eRenderPipeline {
			SINGLEPASS,
			MULTIPASS
		};

		enum eTextureType {
			COMPLETE,
			NORMAL,
			OCCLUSION,
			EMISSIVE,
			METALNESS,
			ROUGHNESS
		};

		enum ePipeline {
			FORWARD,
			DEFERRED
		};
		
		enum eSSAOType {
			SSAO,
			SSAO_plus
		};


		// Save only the visible nodes sorted by distance to the camera
		std::vector<RenderCall> render_calls;
		// Save all lights in the scene
		std::vector<LightEntity*> lights;

		int max_lights;
		int num_lights_shadows;

		// FBOs
		FBO* fbo;
		FBO* gbuffers_fbo;
		FBO* illumination_fbo;
		FBO* screen_fbo;
		FBO* ssao_fbo;
		Texture* shadowmap;
		Texture* screen_texture;
		int width_shadowmap; 
		int height_shadowmap;

		// Imgui debug parameters
		bool show_shadowmap;
		int debug_shadowmap;
		int debug_texture;
		bool show_buffers;
		bool show_ssao;
		bool pbr;

		ePipeline pipeline;
		int typeOfRender;

		// color correction
		bool gamma;
		bool tonemapping;
		float lumwhite2;
		float averagelum;
		float scale;

		// ambient occlusion
		std::vector<Vector3> random_points_sph;
		std::vector<Vector3> random_points_hemi;
		eSSAOType SSAOType;
		bool ssao;

		sProbe probe;
			
		// meshes
		Mesh* sphere;
		Mesh* quad;

		Renderer();

		// -- Rendercalls manager functions--

		void createRenderCalls(GTR::Scene* scene, Camera* camera);
		void addRenderCall_node(GTR::Scene* scene, Camera* camera, Node* node, Matrix44 curr_model, Matrix44 parent_model);
		void sortRenderCalls();
		// operator used to sort rendercalls vector
		static bool compare_distances(const RenderCall rc1, const RenderCall rc2) { return (rc1.distance_to_camera < rc2.distance_to_camera); }

		// -- Shadowmap functions --

		Vector4 assignMapPiece(int width, int height, int index, int num_elements);
		Vector4 assignMapPiece_shader(int width, int height, int index, int num_elements);
		void generateShadowmap(LightEntity* light, int index);
		void showShadowmap(LightEntity* light);

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
		void renderTransparentMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		//Texture* screenToTexture();
		void renderMeshWithMaterialToGBuffers(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		// to render flat objects for generating the shadowmaps
		void renderFlatMesh(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		void renderForward(Camera* camera);
		void renderDeferred(Camera* camera);

		// -- Deferred functions

		void generateGBuffers(Scene* scene, Camera* camera);
		void generateSSAO(Camera* camera, Matrix44 inv_vp, int width, int height);
		void applyIllumination_deferred(Scene* scene, Camera* camera, Matrix44 inv_vp, int width, int height);
		void applyColorCorrection();

		// -- Upload to shader functions --

		void uploadLightToShader(GTR::LightEntity* light, Shader* shader, Vector3 ambient_light, int shadow_i);
		void setTextures(GTR::Material* material, Shader* shader);
		void setSinglepass_parameters(GTR::Material* material, Shader* shader, Mesh* mesh);
		void setMultipassParameters(GTR::Material* material, Shader* shader, Mesh* mesh);

		// -- Debug functions --

		void showGBuffers(int width, int height, Camera* camera);
		void setLightsVisible();
		void setLightsInvisible();
		void renderInMenu();
		
		void renderProbe(Vector3 pos, float size, float* coeffs);
	};

	Texture* CubemapFromHDRE(const char* filename);
	std::vector<Vector3> generateSpherePoints(int num, float radius, bool hemi);
};