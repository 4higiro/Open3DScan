#ifndef _OPEN_3D_SCAN_CALCULATE_3D
#define _OPEN_3D_SCAN_CALCULATE_3D

#include "linear_math.hpp"
#include "image_processing.hpp"
#include <vector>

namespace O3DS
{
	class img_t;

	namespace proc_3d
	{
		enum class camera_selection_t { A, B };

		struct camera_data_t final
		{
			const img_t& source;
			const proc_2d::rectangle_detector_t detector;
			math::vec2 pixel_size;
		};

		struct calibration_data_t final
		{
			double width = 1.0;
			double height = 1.0;
			double depth = 1.0;
		};

		struct cameras_proc_data_t final
		{
			math::vec3 cama_euler, camb_euler;
			math::vec3 cama_position, camb_position;
			double cama_h = 0.0, camb_h = 0.0;
		};

		class cameras_pair_t final
		{
		private:
			double _cama_h = 1.0, _camb_h = 1.0;
			math::size2 _cama_size, _camb_size;
			math::vec2 _cama_pixel_size, _camb_pixel_size;
			bool _calibration_status = false;
		
			static void _calc_sides(const std::array<math::vec2, 4>& points, double& out_a, double& out_b) noexcept;
			static bool _calc_calibration_params(const img_t& source, double width, double height, double depth,
				const proc_2d::rectangle_detector_t& detector, const math::vec2 pixel_size, double& out_h) noexcept;
		public:
			cameras_pair_t(const camera_data_t& cam_a, const camera_data_t& cam_b, const calibration_data_t& data) noexcept;
			cameras_pair_t(double cama_h, math::size2 cama_size, double camb_h, math::size2 camb_size) noexcept;

			bool get_calibration_status() const noexcept;

			double get_camera_h(camera_selection_t pick) const noexcept;
			
			void get_camera_rotation(proc_2d::rectangle_detector_t::rect_t source, math::vec3& target,
				camera_selection_t pick = camera_selection_t::A) const noexcept;

			void get_cameras_position(std::vector<proc_2d::rectangle_detector_t::rect_t> rects_a,
				std::vector<proc_2d::rectangle_detector_t::rect_t> rects_b,
				math::vec3 cama_euler, math::vec3 camb_euler, math::vec3& cama_position, math::vec3& camb_position) noexcept;
		};

		class processor_t final
		{
		private:
			cameras_proc_data_t _cameras;
			proc_2d::line_detector_t _detector;
			double _dist_threshold = math::positive_inf;
		public:
			processor_t(const cameras_proc_data_t& cameras, const proc_2d::line_detector_t& detector,
				double dist_threshold = math::positive_inf) noexcept;

			cameras_proc_data_t get_cameras_proc_data() const noexcept;
			void set_cameras_proc_data(const cameras_proc_data_t& cameras) noexcept;

			proc_2d::line_detector_t get_line_detector() const noexcept;
			void set_line_detector(const proc_2d::line_detector_t& detector) noexcept;

			// safe exist points
			void expand_point_cloud(const img_t& a, const img_t& b, std::vector<math::vec3> point_cloud) const noexcept;
		};
	}
}

#endif;
