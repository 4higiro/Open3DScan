#ifndef _OPEN_3D_SCAN_IMAGE_PROCESSING
#define _OPEN_3D_SCAN_IMAGE_PROCESSING

#include "numerics.hpp"
#include "linear_math.hpp"
#include <vector>
#include <array>
#include <list>

namespace O3DS
{
	namespace compute { extern bool use_gpu; }

	class img_t;

	namespace proc_2d
	{
		class filter1D_t final
		{
		private:
			struct {
				float _T = 0;
				size_t _ra = 0;
				size_t _rb = 0;
			} _coord_system;

			std::vector<float> _arr;
			size_t _arr_size = 0;

		public:
			template <typename function_t>
			filter1D_t(function_t f, size_t ra, size_t rb = math::inf, float T = 0) noexcept
				: _coord_system({ T, ra, rb == math::inf ? ra : rb })
			{
				if (_coord_system._ra == math::inf || _coord_system._rb == math::inf) return;
				_arr.resize(_coord_system._ra + _coord_system._rb + 1);
				_arr_size = _arr.size();
				for (int i = 0; i < _arr_size; ++i)
				{
					double t = -_coord_system._T - _coord_system._ra + i;
					_arr[i] = f(t);
				}
			}

			filter1D_t(std::vector<float>&& arr) noexcept;

			float& operator()(int x) noexcept;
			float operator()(int x) const noexcept;

			size_t size() const noexcept;
			int leftmost() const noexcept;
			int rightmost() const noexcept;
			float bias() const noexcept;
			const float* data() const noexcept;
		};

		class filter2D_t final
		{
		private:
			struct {
				math::vec2f _T;
				math::size2 _rx;
				math::size2 _ry;
			} _coord_system;

			std::vector<float> _arr;
			size_t _width, _height;

		public:
			template <typename function_t>
			filter2D_t(function_t f, const math::size2& rx, const math::size2& ry = { math::inf, math::inf },
				const math::vec2f& T = { 0.0, 0.0 }) noexcept : _coord_system({T, rx, ry["x"] == math::inf || ry["y"] == math::inf ? rx : ry})
			{
				if (_coord_system._rx["x"] == math::inf || _coord_system._rx["y"] == math::inf
					|| _coord_system._ry["x"] == math::inf || _coord_system._ry["y"] == math::inf)
					return;
				_width = _coord_system._rx["x"] + _coord_system._rx["y"] + 1;
				_height = _coord_system._ry["x"] + _coord_system._ry["y"] + 1;
				_arr.resize(_width * _height);
				for (int i = 0; i < _height; ++i)
				{
					for (int j = 0; j < _width; ++j)
					{
						double x = -_coord_system._T["x"] - _coord_system._rx["x"] + j;
						double y = -_coord_system._T["y"] - _coord_system._ry["x"] + i;
						_arr[j + i * _width] = f(x, y);
					}
				}
			}

			filter2D_t(std::vector<float>&& arr, size_t w, size_t h) noexcept;

			float& operator()(int x, int y) noexcept;
			float operator()(int x, int y) const noexcept;

			size_t width() const noexcept;
			size_t height() const noexcept;
			int leftmost() const noexcept;
			int rightmost() const noexcept;
			int topmost() const noexcept;
			int lowermost() const noexcept;
			const float* data() const noexcept;
		};

		struct img_row_t final { const img_t& data; size_t row_number; };

		float convolution1D(const img_row_t& img, const filter1D_t& filter, size_t x) noexcept;
		float convolution2D(const img_t& img, const filter2D_t& filter, const math::size2& x) noexcept;
		math::vec2f sobel_operator(const img_t& img, const math::size2& x) noexcept;
		img_t canny_edge_detection(const img_t& img, double sigma, float threshold) noexcept;

		enum line_mode_t { VERTICAL, HORIZONTAL, UNDEFINE = VERTICAL };

		class refiner_interface_t
		{
		public:
			virtual void init(double sigma, double min_threshold) noexcept = 0;
			virtual bool refine(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept = 0;
		};

		class avgweight_refiner_t : public refiner_interface_t
		{
		private:
			int _r = 1;

		public:
			void init(double sigma, double min_threshold) noexcept override;
			bool refine(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept override;

			avgweight_refiner_t();
			avgweight_refiner_t(double sigma) noexcept;
		};

		class stager_1D_refiner_t : public refiner_interface_t
		{
		private:
			double _sigma = 1.5;
			double _max_threshold = 2.0;
			int _filter_semi_len = 4;

			double _refine_iteration(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept;
		public:
			static constexpr int iteration_count = 10;
			static constexpr double eps = 1e-4;

			void init(double sigma, double max_threshold) noexcept override;
			bool refine(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept override;

			stager_1D_refiner_t();
			stager_1D_refiner_t(double sigma, double max_threshold = 2.0) noexcept;
		};

		class stager_2D_refiner_t : public refiner_interface_t
		{
		private:
			double _sigma = 1.5;
			double _min_threshold = 2.0;
			int _filter_semi_len = 4;

			double _refine_iteration(const img_t& img, math::vec2& pix) const noexcept;
		public:
			static constexpr int iteration_count = 10;
			static constexpr double eps = 1e-4;

			void init(double sigma, double min_threshold) noexcept override;
			bool refine(const img_t& img, math::vec2& pix, line_mode_t) const noexcept override;

			stager_2D_refiner_t();
			stager_2D_refiner_t(double sigma, double min_threshold = 2.0) noexcept;
		};

		class line_detector_t final
		{
		private:
			line_mode_t _mode = UNDEFINE;
			std::array<math::size2, 2> _border;
			refiner_interface_t* _refiner = nullptr;
			uint8_t _intensity_threshold = 100;

			uint8_t _find_max_value(const img_row_t& img, size_t a, size_t b, size_t& max_index) const noexcept;
		public:
			line_detector_t(refiner_interface_t* refiner, uint8_t intensity_threshold = 100,
				line_mode_t mode = line_mode_t::UNDEFINE, math::size2 left_top_point = { math::inf, math::inf },
				math::size2 right_buttom_point = { math::inf, math::inf }) noexcept;

			line_mode_t get_mode() const noexcept;
			void set_mode(line_mode_t mode) noexcept;

			std::array<math::size2, 2> get_border() const noexcept;
			void set_border(const std::array<math::size2, 2>& border) noexcept;

			uint8_t get_intensity_threshold() const noexcept;
			void set_intensity_threshold(uint8_t intensity_threshold) noexcept;

			void set_refiner(refiner_interface_t* refiner) noexcept;

			void get_line(const img_t& img, std::list<math::vec2>& points) const noexcept; // safe existest points
			math::vec2 get_point(const img_t& img, math::vec2 pixel) const noexcept;
		};

		class rectangle_detector_t final
		{
		private:
			float _canny_threshold = 75.0f;
			float _match_threshold = 0.5f;
			double _equal_threshold = 10.0;
			size_t _angle_size = 15;
			size_t _nearest_radius = 5;
			size_t _line_semi_width;

			void _gpu_find_angles(const img_t& img, uint8_t*& result) const noexcept;
			void _cpu_find_angles(const O3DS::img_t& temp, uint8_t*& angle_mask) const noexcept;

		public:
			struct rect_t final
			{
				std::array<math::vec2, 4> verts;

				double square() const noexcept;
				static bool equal(const rect_t& a, const rect_t& b, double threshold) noexcept;
			};

			rectangle_detector_t(size_t angle_size = 15, float canny_threshold = 75.0f, size_t line_semi_width = 7, size_t nearest_radius = 5,
				float match_threshold = 0.5f, double equal_threshold = 10.0) noexcept;

			float get_canny_threshold() const noexcept;
			void set_canny_threshold(float threshold) noexcept;

			size_t get_angle_size() const noexcept;
			void set_angle_size(size_t angle_size) noexcept;

			float get_match_threshold() const noexcept;
			void set_match_threshold(float threshold) noexcept;

			float get_equal_threshold() const noexcept;
			void set_equal_threshold(float threshold) noexcept;

			float get_nearest_radius() const noexcept;
			void set_nearest_radius(size_t radius) noexcept;
			
			float get_line_semi_width() const noexcept;
			void set_line_semi_width(size_t semi_width) noexcept;

			void get_rectangles(const img_t& img, std::vector<rect_t>& rects) const noexcept;
		};
	}
}

#endif
