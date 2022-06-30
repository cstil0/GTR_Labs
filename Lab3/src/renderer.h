#pragma once
#include "prefab.h"
#include "shader.h"
#include "sphericalharmonics.h"
#include "mesh.h"
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

	struct sReflectionProbe {
		Vector3 pos;
		float size = 10;
		Texture* cubemap = NULL;
	};


	struct sIrrHeader {
		Vector3 start;
		Vector3 end;
		Vector3 delta;
		Vector3 dims;
		int num_probes;
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
		std::vector<DecalEntity*> decals;
		LightEntity* direct_light;

		int max_lights;
		int num_lights_shadows;

		// FBOs
		FBO* fbo;
		FBO* gbuffers_fbo;
		FBO* illumination_fbo;
		FBO* irr_fbo;
		FBO* screen_fbo;
		FBO* ssao_fbo;
		FBO* decals_fbo;
		FBO* volumetric_fbo;
		FBO* linear_depth_fbo;
		Texture* shadowmap;
		Texture* screen_texture;
		Texture* probes_texture;
		Texture* postFX_textureA;
		Texture* postFX_textureB;
		Texture* postFX_textureC;
		Texture* postFX_textureD;
		Texture* randomTexture;


		// debuging use -- borrar
		Texture* irradiance_texture;
		int width_shadowmap; 
		int height_shadowmap;

		// Imgui debug parameters
		bool show_shadowmap;
		int debug_shadowmap;
		int debug_texture;
		bool show_buffers;
		bool show_ssao;
		bool pbr;
		bool global_metrough;

		float global_metalness;
		float global_roughness;

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
		bool show_probes;
		bool show_probes_texture;
		bool shadow_flag;

		// irradiance
		std::vector<sProbe> irradiance_probes;
		Vector3 start_irr;
		Vector3 end_irr;
		Vector3 dim_irr;
		Vector3 delta_irr;
		bool irradiance;
		bool show_irradiance;

		// meshes
		Mesh* sphere;
		Mesh* quad;
		Mesh cube;
		Mesh* cone;

		// skybox
		//Texture* skybox;

		// reflections
		std::vector<sReflectionProbe*> reflection_probes;
		FBO* reflection_fbo;
		FBO* reflection_probe_fbo;
		bool is_rendering_reflections;
		bool is_rendering_planar_reflections;
		bool planar_reflection;
		bool render_reflection_probes;
		bool scene_reflection;
		//ReflectionProbeEntity* reflection_probe;

		// volumetric
		bool volumetric;

		// PostFX
		bool postFX;
		bool saturation;
		bool lens_distortion;
		bool contrast;
		bool simple_glow;
		bool perfect_glow;
		bool depth_field;
		bool antialiasing;
		bool grain;

		float saturation_intensity;
		float vigneting_intensity;
		float contrast_intensity;
		float simglow_blur_factor;
		float simglow_mix_factor;
		float simglow_threshold;
		int perfglow_iterations;
		float apperture;
		float focal_length;
		float focal_range;
		bool show_depth_field;
		float grainIntensity;

		Matrix44 vp_matrix_last;

		Renderer();

		// -- Rendercalls manager functions --

		void createRenderCalls(Camera* camera);
		void addRenderCall_node(Camera* camera, Node* node, Matrix44 curr_model, Matrix44 parent_model);
		void sortRenderCalls();
		// operator used to sort rendercalls vector
		static bool compare_distances(const RenderCall rc1, const RenderCall rc2) { return (rc1.distance_to_camera < rc2.distance_to_camera); }

		// -- Shadowmap functions --

		Vector4 assignMapPiece(int width, int height, int index, int num_elements);
		Vector4 assignMapPiece_shader(int width, int height, int index, int num_elements);
		void generateShadowmap(LightEntity* light, int index);
		void generateIrradianceProbes();
		void showShadowmap(LightEntity* light);

		void lineralizeDepth(Texture* depth, Vector2 camera_nearfar);

		// -- Render functions --
		
		//renders several elements of the scene
		void renderScene(GTR::Scene* scene, Camera* camera);
		void renderSceneWithReflection(Camera* camera);
		// to render the scene using rendercalls vector
		void renderScene_RenderCalls(Camera* camera, FBO* fboToRender);
		//to render a whole prefab (with all its nodes)
		void renderPrefab(const Matrix44& model, GTR::Prefab* prefab, Camera* camera);
		//to render one node from the prefab and its children
		void renderNode(const Matrix44& model, GTR::Node* node, Camera* camera);
		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera, int i);
		void renderTransparentMaterial(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		//Texture* screenToTexture();
		void renderMeshWithMaterialToGBuffers(Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		// to render flat objects for generating the shadowmaps
		void renderFlatMesh(const Matrix44 model, Mesh* mesh, GTR::Material* material, Camera* camera);
		void renderForward(Camera* camera, FBO* fboToRender);
		void renderDeferred(Camera* camera);

		void renderSkybox(Camera* camera, Texture* skybox);

		// -- Deferred functions

		void generateGBuffers(Camera* camera);
		void generateSSAO(Camera* camera, Matrix44 inv_vp, int width, int height);
		void applyIllumination_deferred(Camera* camera, Matrix44 inv_vp, int width, int height);
		void applyColorCorrection(Texture* final_texture);

		// -- Upload to shader functions --

		void uploadLightToShader(GTR::LightEntity* light, Shader* shader, Vector3 ambient_light, int shadow_i = NULL);
		void setTextures(GTR::Material* material, Shader* shader);
		void setSinglepass_parameters(GTR::Material* material, Shader* shader, Mesh* mesh);
		void setMultipassParameters(GTR::Material* material, Shader* shader, Mesh* mesh);

		// -- Debug functions --

		void showGBuffers(int width, int height, Camera* camera);
		void setLightsVisible();
		void setLightsInvisible();
		void renderInMenu();
		
		void renderIrradianceProbe(Vector3 pos, float size, float* coeffs);
		void captureIrradianceProbe(sProbe& probe);
		void saveIrradianceProbesToDisk();
		bool loadIrradianceProbesFromDisk();
		void computeIrradianceDeferred(Matrix44 inv_vp);
		void computeIrradianceForward(Mesh* mesh, Matrix44 model, Material* material, Camera* camera, int i);
		
		void generateReflectionProbes();
		void generateReflectionProbesMesh();
		void renderReflectionProbes(Camera* camera);
		void captureReflectionProbe(Texture* tex, Vector3 pos);
		Texture* applyAntialiasing(Shader* fxshader, Texture* draw_texture, Texture* read_texture);
		Texture* applySaturation(Shader* fxshader, Texture* draw_texture, Texture* read_texture);
		Texture* applyLensDistortion(Shader* fxshader, Texture* draw_texture, Texture* read_texture);
		Texture* applyContrast(Shader* fxshader, Texture* draw_texture, Texture* read_texture);
		Texture* applySimpleGlow(Shader* fxshader, Texture* draw_texture1, Texture* draw_texture2, Texture* draw_texture3, Texture* draw_texture4, Texture* read_texture);
		Texture* applyPerfectGlow(Shader* fxshader, Texture* draw_texture1, Texture* draw_texture2, Texture* draw_texture3, Texture* draw_texture4, Texture* read_texture);
		Texture* applyBlurring(Shader* fxshader, Texture* draw_texture1, Texture* draw_texture2, Texture* read_texture, int iterations);
		Texture* applyDepthField(Shader* fxshader, Camera* camera, Texture* draw_texture1, Texture* draw_texture2, Texture* draw_texture3, Texture* draw_texture4, Texture* read_texture, Texture* depth_texture);
		Texture* applyGrain(Shader* fxshader, Camera* camera, Texture* draw_texture, Texture* read_texture);
		void saveReflectionProbesToDisk();
		bool loadReflectionProbesFromDisk();

		Texture* applyFX(Camera* camera, Texture* color_texture, Texture* depth_texture);

		// FUERA
		void reloadRenderer();
};

	std::vector<Vector3> generateSpherePoints(int num, float radius, bool hemi);
};