#include "image.hpp"
#include "image_processing.hpp"
#include <algorithm>

#define BOOST_CHRONO_NO_LIB
#include <boost/compute.hpp>

namespace O3DS { namespace compute {
	extern bool use_gpu;
}}

static boost::compute::device gpu_device = boost::compute::system::default_device();

O3DS::proc_2d::filter1D_t::filter1D_t(std::vector<float>&& arr) noexcept
	: _arr(arr), _arr_size(arr.size()), _coord_system({ 0.0, arr.size() / 2, arr.size() / 2 }) {}

float& O3DS::proc_2d::filter1D_t::operator()(int x) noexcept
{
	return _arr[x + _coord_system._ra];
}

float O3DS::proc_2d::filter1D_t::operator()(int x) const noexcept
{
	return _arr[x + _coord_system._ra];
}

size_t O3DS::proc_2d::filter1D_t::size() const noexcept
{
	return _arr_size;
}

int O3DS::proc_2d::filter1D_t::leftmost() const noexcept
{
	return -static_cast<int>(_coord_system._ra);
}

int O3DS::proc_2d::filter1D_t::rightmost() const noexcept
{
	return static_cast<int>(_coord_system._rb);
}

float O3DS::proc_2d::filter1D_t::bias() const noexcept
{
	return _coord_system._T;
}

const float* O3DS::proc_2d::filter1D_t::data() const noexcept
{
	return _arr.data();
}

O3DS::proc_2d::filter2D_t::filter2D_t(std::vector<float>&& arr, size_t w, size_t h) noexcept
	: _arr(arr), _width(w), _height(h), _coord_system({ { 0.0, 0.0 }, { w / 2, w / 2 }, { h / 2, h / 2 } }) {}

float& O3DS::proc_2d::filter2D_t::operator()(int x, int y) noexcept
{
	return _arr[_coord_system._rx["x"] + x + (_coord_system._ry["x"] + y) * _width];
}

float O3DS::proc_2d::filter2D_t::operator()(int x, int y) const noexcept
{
	return _arr[_coord_system._rx["x"] + x + (_coord_system._ry["x"] + y) * _width];
}

size_t O3DS::proc_2d::filter2D_t::width() const noexcept
{
	return _width;
}

size_t O3DS::proc_2d::filter2D_t::height() const noexcept
{
	return _height;
}

int O3DS::proc_2d::filter2D_t::leftmost() const noexcept
{
	return -static_cast<int>(_coord_system._rx["x"]);
}

int O3DS::proc_2d::filter2D_t::rightmost() const noexcept
{
	return static_cast<int>(_coord_system._rx["y"]);
}

int O3DS::proc_2d::filter2D_t::topmost() const noexcept
{
	return static_cast<int>(_coord_system._ry["y"]);
}

int O3DS::proc_2d::filter2D_t::lowermost() const noexcept
{
	return -static_cast<int>(_coord_system._ry["x"]);
}

const float* O3DS::proc_2d::filter2D_t::data() const noexcept
{
	return _arr.data();
}

float O3DS::proc_2d::convolution1D(const img_row_t& img, const filter1D_t& filter, size_t arg) noexcept
{
	float integral = 0.0;
	for (int x = filter.leftmost(); x <= filter.rightmost(); ++x)
		integral += img.data(arg - x, img.row_number) * filter(x);
	return integral;
}

float O3DS::proc_2d::convolution2D(const img_t& img, const filter2D_t& filter, const math::size2& arg) noexcept
{
	float integral = 0.0;
	for (int y = filter.lowermost(); y <= filter.topmost(); ++y)
	{
		for (int x = filter.leftmost(); x <= filter.rightmost(); ++x)
			integral += img(arg["x"] - x, arg["y"] - y) * filter(x, y);
	}
	return integral;
}

O3DS::math::vec2f O3DS::proc_2d::sobel_operator(const img_t& img, const math::size2& arg) noexcept
{
	std::vector<float> mgx = {
		1.0, 0.0, -1.0,
		2.0, 0.0, -2.0,
		1.0, 0.0, -1.0
	};
	std::vector<float> mgy = {
		1.0, 2.0, 1.0,
		0.0, 0.0, 0.0,
		-1.0, -2.0, -1.0
	};
	filter2D_t filter_x(std::move(mgx), 3, 3);
	filter2D_t filter_y(std::move(mgy), 3, 3);
	float x = convolution2D(img, filter_x, arg);
	float y = convolution2D(img, filter_y, arg);
	return { x, y };
}

size_t get_range_size(size_t size, size_t work_group_size) noexcept
{
	size_t n = size / work_group_size + 1;
	size_t mod = n * work_group_size - size;
	return size + mod;
}

void gpu_canny_edge_detection(O3DS::img_t& img, const O3DS::proc_2d::filter2D_t& gauss_filter, float threshold)
{
	static constexpr const char* convolution2D_kernel_source = R"(
		__kernel void conv(__global uchar* in, __global uchar* out, __global float* filter, 
			__global uint* pfilter_semi_len, __global uint* pw, __global uint* ph)
		{
			int filter_semi_len = pfilter_semi_len[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);
			
			if (index < w * h)
			{
				int filter_len = filter_semi_len * 2 + 1;
				float integral = 0.0;
				for (int i = -filter_semi_len; i <= filter_semi_len; ++i)
				{
					for (int j = -filter_semi_len; j <= filter_semi_len; ++j)
					{
						if (j0 - j < 0 || j0 - j >= w || i0 - i < 0 || i0 - i >= h)
							continue;
						uint zi = j0 - j + (i0 - i) * w;
						uint fy = i + filter_semi_len;
						uint fx = j + filter_semi_len;
						integral += convert_float(in[zi]) * filter[fx + fy * filter_len];
					}
				}

				out[j0 + i0 * w] = convert_uchar(integral);
			}
		}
	)";

	static constexpr const char* fconvolution2D_kernel_source = R"(
		__kernel void fconv(__global uchar* in, __global float* out, __global float* filter, 
			__global uint* pfilter_semi_len, __global uint* pw, __global uint* ph)
		{
			int filter_semi_len = pfilter_semi_len[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				int filter_len = filter_semi_len * 2 + 1;
				float integral = 0.0;
				for (int i = -filter_semi_len; i <= filter_semi_len; ++i)
				{
					for (int j = -filter_semi_len; j <= filter_semi_len; ++j)
					{
						if (j0 - j < 0 || j0 - j >= w || i0 - i < 0 || i0 - i >= h)
							continue;
						uint zi = j0 - j + (i0 - i) * w;
						uint fy = i + filter_semi_len;
						uint fx = j + filter_semi_len;
						integral += convert_float(in[zi]) * filter[fx + fy * filter_len];
					}
				}

				out[j0 + i0 * w] = integral;
			}
		}
	)";

	static constexpr const char* detection_kernel_source = R"(
		__kernel void detect(__global float* px, __global float* py, __global uchar* out,
			__global uint* pw, __global uint* ph, __global float* pthreshold)
		{
			int w = pw[0];
			int h = ph[0];
			float threshold = pthreshold[0]; 

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				float vx = px[j0 + i0 * w];
				float vy = py[j0 + i0 * w];
				vx = vx / sqrt(vx * vx + vy * vy);
				vy = vy / sqrt(vx * vx + vy * vy);

				int ix = 0, iy = 0;
				if (vx > 0.5f)	ix = 1;
				if (vx < -0.5f) ix = -1;
				if (vy > 0.5f)	iy = 1;
				if (vy < -0.5f) iy = -1;

				if (j0 + ix >= 0 && j0 - ix >= 0 && j0 + ix < w && j0 - ix < w
					&& i0 + iy >= 0 && i0 - iy >= 0 && i0 + iy < h && i0 - iy < h)
				{
					float v1x = px[j0 + ix + (i0 + iy) * w];
					float v1y = py[j0 + ix + (i0 + iy) * w];
					float v2x = px[j0 - ix + (i0 - iy) * w];
					float v2y = py[j0 - ix + (i0 - iy) * w];
					float v0x = px[j0 + i0 * w];
					float v0y = py[j0 + i0 * w];
					float l1 = v1x * v1x + v1y * v1y;
					float l2 = v2x * v2x + v2y * v2y;
					float l0 = v0x * v0x + v0y * v0y;
					if (l0 > l1 && l0 > l2 && l0 > threshold)
						out[j0 + i0 * w] = 255;
					else
						out[j0 + i0 * w] = 0;
				}
				else
					out[j0 + i0 * w] = 0;
			}
		}
	)";

	boost::compute::context gpu_context(gpu_device);
	boost::compute::command_queue gpu_command_queue(gpu_context, gpu_device);

	const uint8_t* img_data = img.get_format()->bitmap;
	uint8_t* output_data = const_cast<uint8_t*>(img_data);
	cl_uint w = img.width();
	cl_uint h = img.height();
	boost::compute::program conv_program = boost::compute::program::create_with_source(convolution2D_kernel_source, gpu_context);
	conv_program.build();
	boost::compute::kernel conv_kernel(conv_program, "conv");
	boost::compute::program fconv_program = boost::compute::program::create_with_source(fconvolution2D_kernel_source, gpu_context);
	fconv_program.build();
	boost::compute::kernel fconv_kernel(fconv_program, "fconv");
	boost::compute::program detect_program = boost::compute::program::create_with_source(detection_kernel_source, gpu_context);
	detect_program.build();
	boost::compute::kernel detect_kernel(detect_program, "detect");
	boost::compute::buffer input(gpu_context, w * h);
	boost::compute::buffer blur(gpu_context, w * h);
	boost::compute::buffer fslb(gpu_context, sizeof(cl_uint));
	boost::compute::buffer wb(gpu_context, sizeof(cl_uint));
	boost::compute::buffer hb(gpu_context, sizeof(cl_uint));
	gpu_command_queue.enqueue_write_buffer(wb, 0, sizeof(cl_uint), &w);
	gpu_command_queue.enqueue_write_buffer(hb, 0, sizeof(cl_uint), &h);

	// Gaussing blur
	{
		const float* filter_data = gauss_filter.data();
		cl_uint filter_semi_len = gauss_filter.width() / 2;
		size_t filter_buff_size = gauss_filter.width() * gauss_filter.height() * sizeof(float);
		boost::compute::buffer filter(gpu_context, filter_buff_size);
		gpu_command_queue.enqueue_write_buffer(input, 0, w * h, img_data);
		gpu_command_queue.enqueue_write_buffer(filter, 0, filter_buff_size, filter_data);
		gpu_command_queue.enqueue_write_buffer(fslb, 0, sizeof(cl_uint), &filter_semi_len);
		conv_kernel.set_arg(0, input);
		conv_kernel.set_arg(1, blur);
		conv_kernel.set_arg(2, filter);
		conv_kernel.set_arg(3, fslb);
		conv_kernel.set_arg(4, wb);
		conv_kernel.set_arg(5, hb);
		gpu_command_queue.enqueue_1d_range_kernel(conv_kernel, 0, get_range_size(w * h, 128), 128);
	}

	boost::compute::buffer der_x(gpu_context, w * h * sizeof(cl_float));
	boost::compute::buffer der_y(gpu_context, w * h * sizeof(cl_float));

	// Sobel operations
	{
		std::vector<float> mgx = {
		1.0, 0.0, -1.0,
		2.0, 0.0, -2.0,
		1.0, 0.0, -1.0
		};
		const float* filter_x_data = mgx.data();
		cl_uint filter_semi_len = 1;
		size_t filter_buff_size = 9 * sizeof(float);
		boost::compute::buffer filter(gpu_context, filter_buff_size);
		gpu_command_queue.enqueue_write_buffer(filter, 0, filter_buff_size, filter_x_data);
		gpu_command_queue.enqueue_write_buffer(fslb, 0, sizeof(cl_uint), &filter_semi_len);
		fconv_kernel.set_arg(0, blur);
		fconv_kernel.set_arg(1, der_x);
		fconv_kernel.set_arg(2, filter);
		fconv_kernel.set_arg(3, fslb);
		fconv_kernel.set_arg(4, wb);
		fconv_kernel.set_arg(5, hb);
		gpu_command_queue.enqueue_1d_range_kernel(fconv_kernel, 0, get_range_size(w * h, 128), 128);
		std::vector<float> mgy = {
			1.0, 2.0, 1.0,
			0.0, 0.0, 0.0,
			-1.0, -2.0, -1.0
		};
		const float* filter_y_data = mgy.data();
		gpu_command_queue.enqueue_write_buffer(filter, 0, filter_buff_size, filter_y_data);
		fconv_kernel.set_arg(0, blur);
		fconv_kernel.set_arg(1, der_y);
		fconv_kernel.set_arg(2, filter);
		fconv_kernel.set_arg(3, fslb);
		fconv_kernel.set_arg(4, wb);
		fconv_kernel.set_arg(5, hb);
		gpu_command_queue.enqueue_1d_range_kernel(fconv_kernel, 0, get_range_size(w * h, 128), 128);
	}

	boost::compute::buffer& output = input;

	// Edge detection
	{
		boost::compute::buffer thrb(gpu_context, sizeof(cl_float));
		gpu_command_queue.enqueue_write_buffer(thrb, 0, sizeof(cl_float), &threshold);
		detect_kernel.set_arg(0, der_x);
		detect_kernel.set_arg(1, der_y);
		detect_kernel.set_arg(2, output);
		detect_kernel.set_arg(3, wb);
		detect_kernel.set_arg(4, hb);
		detect_kernel.set_arg(5, thrb);
		gpu_command_queue.enqueue_1d_range_kernel(detect_kernel, 0, get_range_size(w * h, 128), 128);
	}

	gpu_command_queue.enqueue_read_buffer(output, 0, w * h, output_data);
}

void cpu_canny_edge_detection(O3DS::img_t& result, size_t filter_semi_len,
	const O3DS::proc_2d::filter2D_t& gauss_filter, float threshold) noexcept
{
	size_t w = result.width();
	size_t h = result.height();
	for (int i = filter_semi_len; i < h - filter_semi_len; ++i)
	{
		for (int j = filter_semi_len; j < w - filter_semi_len; ++j)
			result(j, i) = O3DS::proc_2d::convolution2D(result, gauss_filter, { j, i });
	}
	std::vector<std::pair<O3DS::math::vec2f, float>> gradient_img(w * h, std::make_pair(O3DS::math::vec2f(), 0));
	for (int i = filter_semi_len; i < h - filter_semi_len; ++i)
	{
		for (int j = filter_semi_len; j < w - filter_semi_len; ++j)
		{
			O3DS::math::vec2f gradient = O3DS::proc_2d::sobel_operator(result, { j, i });
			static const float ang45 = sqrt(2.0) / 2.0;
			static const std::vector<O3DS::math::vec2f> responses = { { 1.0f, 0.0f }, { ang45, ang45 },
																{ 0.0f, 1.0f }, { -ang45, ang45 },
																{ -1.0f, 0.0f }, { -ang45, -ang45 },
																{ 0.0f, -1.0f }, { ang45, -ang45 } };
			gradient_img[j + i * w] = std::make_pair(*std::max_element(responses.begin(), responses.end(),
				[&gradient](const O3DS::math::vec2f& a, const O3DS::math::vec2f& b) 
				{ return O3DS::math::dot(gradient, a) < O3DS::math::dot(gradient, b); }), O3DS::math::length(gradient));
		}
	}
	for (int i = 0; i < filter_semi_len + 1; ++i)
	{
		for (int j = 0; j < w; ++j)
			result(j, i) = 0;
	}
	for (int i = filter_semi_len + 1; i < h - filter_semi_len - 1; ++i)
	{
		for (int j = 0; j < filter_semi_len + 1; ++j)
			result(j, i) = 0;
		for (int j = filter_semi_len + 1; j < w - filter_semi_len - 1; ++j)
		{
			O3DS::math::vec2i prev = O3DS::math::vec2i(j, i) - O3DS::math::round(gradient_img[j + i * w].first);
			O3DS::math::vec2i next = O3DS::math::vec2i(j, i) + O3DS::math::round(gradient_img[j + i * w].first);
			float l1 = gradient_img[prev["x"] + prev["y"] * w].second;
			float l2 = gradient_img[next["x"] + next["y"] * w].second;
			float l0 = gradient_img[j + i * w].second;
			if (l0 > l1 && l0 > l2 && l0 > threshold)
				result(j, i) = 0b1111'1111;
			else
				result(j, i) = 0b0000'0000;
		}
		for (int j = w - filter_semi_len - 1; j < w; ++j)
			result(j, i) = 0;
	}
	for (int i = h - filter_semi_len - 1; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
			result(j, i) = 0;
	}
}

O3DS::img_t O3DS::proc_2d::canny_edge_detection(const img_t& img, double sigma, float threshold) noexcept
{
	img_t result = img.get_grayscale();
	threshold *= sqrt(3.0f * 255.0f * 255.0f);
	const size_t filter_semi_len = 3 * sigma;
	filter2D_t gauss_filter([&sigma](double x, double y)
		{ return math::gauss(x, sigma) * math::gauss(y, sigma); }, { filter_semi_len, filter_semi_len });
	if (compute::use_gpu)
		gpu_canny_edge_detection(result, gauss_filter, threshold);
	else
		cpu_canny_edge_detection(result, filter_semi_len, gauss_filter, threshold);

	return result;
}

uint8_t O3DS::proc_2d::line_detector_t::_find_max_value(const img_row_t& img, size_t a, size_t b, size_t& max_index) const noexcept
{
	uint8_t max_value = 0;
	a = a == math::inf ? 0 : a;
	b = b == math::inf ? (_mode == line_mode_t::VERTICAL ? img.data.width() : img.data.height()) : b;
	max_index = a;
	switch (_mode)
	{
	case line_mode_t::VERTICAL:
		for (size_t i = a; i <= b; ++i)
		{
			if (img.data(i, img.row_number) > max_value)
			{
				max_value = img.data(i, img.row_number);
				max_index = i;
			}
		}
		return max_value;
	case line_mode_t::HORIZONTAL:
		for (size_t i = a; i <= b; ++i)
		{
			if (img.data(img.row_number, i) > max_value)
			{
				max_value = img.data(img.row_number, i);
				max_index = i;
			}
		}
		return max_value;
	default:
		return max_value;
	}
}

O3DS::proc_2d::line_detector_t::line_detector_t(refiner_interface_t* refiner, uint8_t intensity_threshold, line_mode_t mode,
	math::size2 left_top_point, math::size2 right_buttom_point) noexcept
	: _refiner(refiner), _mode(mode), _intensity_threshold(intensity_threshold), _border({left_top_point, right_buttom_point}) {}

O3DS::proc_2d::line_mode_t O3DS::proc_2d::line_detector_t::get_mode() const noexcept
{
	return _mode;
}

void O3DS::proc_2d::line_detector_t::set_mode(line_mode_t mode) noexcept
{
	_mode = mode;
}

std::array<O3DS::math::size2, 2> O3DS::proc_2d::line_detector_t::get_border() const noexcept
{
	return _border;
}

void O3DS::proc_2d::line_detector_t::set_border(const std::array<math::size2, 2>& border) noexcept
{
	_border = border;
}

uint8_t O3DS::proc_2d::line_detector_t::get_intensity_threshold() const noexcept
{
	return _intensity_threshold;
}

void O3DS::proc_2d::line_detector_t::set_intensity_threshold(uint8_t intensity_threshold) noexcept
{
	_intensity_threshold = intensity_threshold;
}

void O3DS::proc_2d::line_detector_t::set_refiner(refiner_interface_t* refiner) noexcept
{
	_refiner = refiner;
}

void O3DS::proc_2d::line_detector_t::get_line(const img_t& img, std::list<math::vec2>& points) const noexcept
{
	switch (_mode)
	{
	case line_mode_t::VERTICAL:
	{
#pragma omp parallel for
		for (int i = _border[0]["y"]; i < _border[1]["y"]; ++i)
		{
			size_t max_index = _border[0]["x"];
			size_t row_number = i;
			if (_find_max_value({ img, row_number }, _border[0]["x"], _border[1]["x"], max_index) < _intensity_threshold)
				continue;
			math::vec2 subpix = { max_index, i };
			if (_refiner->refine(img, subpix, _mode))
			{
#pragma omp critical
				{
					points.push_back(subpix);
				}
			}
		}
		break;
	}
	case line_mode_t::HORIZONTAL:
	{
#pragma omp parallel for
		for (int i = _border[0]["x"]; i < _border[1]["x"]; ++i)
		{
			size_t max_index = _border[0]["y"];
			size_t column_number = i;
			if (_find_max_value({ img, column_number }, _border[0]["y"], _border[1]["y"], max_index) < _intensity_threshold)
				continue;
			math::vec2 subpix = { i, max_index };
			if (_refiner->refine(img, subpix, _mode))
			{
#pragma omp critical
				{
					points.push_back(subpix);
				}
			}
		}
		break;
	}
	}
}

O3DS::math::vec2 O3DS::proc_2d::line_detector_t::get_point(const img_t& img, math::vec2 pixel) const noexcept
{
	math::vec2 refine_pix = pixel;
	if (_refiner->refine(img, refine_pix, _mode))
		return refine_pix;
	else
		return pixel;
}

void O3DS::proc_2d::avgweight_refiner_t::init(double sigma, double min_threshold) noexcept
{
	_r = abs(3.0 * sigma);
}

bool O3DS::proc_2d::avgweight_refiner_t::refine(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept
{
	math::vec2i ipix= pix;
	switch (mode)
	{
	case line_mode_t::VERTICAL:
	{
		if (pix["x"] - _r < 0 || pix["x"] + _r >= img.width())
			return 0.0;
		double avg = 0.0;
		for (int i = ipix["x"] - _r; i < ipix["x"]; ++i)
			avg -= static_cast<double>(img(i, ipix["y"])) / 255.0 / static_cast<double>(_r);
		for (int i = ipix["x"] + 1; i <= ipix["x"] + _r; ++i)
			avg += static_cast<double>(img(i, ipix["y"])) / 255.0 / static_cast<double>(_r);
		pix["x"] += avg;
		return abs(avg) < 1.0;
	}
	case line_mode_t::HORIZONTAL:
	{
		if (pix["y"] - _r < 0 || pix["y"] + _r >= img.height())
			return 0.0;
		double avg = 0.0;
		for (int i = ipix["y"] - _r; i < ipix["y"]; ++i)
			avg -= static_cast<double>(img(ipix["x"], i)) / 255.0 / static_cast<double>(_r);
		for (int i = ipix["y"] + 1; i <= ipix["y"] + _r; ++i)
			avg += static_cast<double>(img(ipix["x"], i)) / 255.0 / static_cast<double>(_r);
		pix["y"] += avg;
		return abs(avg) < 1.0;
	}
	default:
		return false;
	}
}

O3DS::proc_2d::avgweight_refiner_t::avgweight_refiner_t() = default;

O3DS::proc_2d::avgweight_refiner_t::avgweight_refiner_t(double sigma) noexcept { init(sigma, 0); }

double O3DS::proc_2d::stager_1D_refiner_t::_refine_iteration(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept
{
	math::size2 ipix = pix;
	if (ipix["x"] - _filter_semi_len < 0 || ipix["x"] + _filter_semi_len >= img.width())
		return math::positive_inf;
	double T = pix["x"] - ipix["x"];
	filter1D_t gt([this](double x){ return math::d_dt_gauss(x, _sigma); }, _filter_semi_len, _filter_semi_len, -T);
	filter1D_t gtt([this](double x) { return math::d2_dt2_gauss(x, _sigma); }, _filter_semi_len, _filter_semi_len, -T);
	double rx = convolution1D({ img, ipix["y"] }, gt, ipix["x"]);
	double rxx = convolution1D({ img, ipix["y"] }, gtt, ipix["x"]);
	if (abs(rxx) < eps) return 0.0;
	double k = -rx / rxx;
	if (abs(k) > 1.0) return math::positive_inf;
	pix["x"] += k;
	return rxx;
}

void O3DS::proc_2d::stager_1D_refiner_t::init(double sigma, double max_threshold) noexcept
{
	_sigma = sigma;
	_max_threshold = max_threshold;
	_filter_semi_len = 3 * _sigma;
}

bool O3DS::proc_2d::stager_1D_refiner_t::refine(const img_t& img, math::vec2& pix, line_mode_t mode) const noexcept
{
	math::vec2i ipix = pix;
	switch (mode)
	{
	case line_mode_t::VERTICAL:
	{
		math::vec2 prev = pix;
		double result = 0.0;
		for (int i = 0; i < iteration_count; ++i)
		{
			result = _refine_iteration(img, pix, mode);
			if (result > _max_threshold)
			{
				pix = prev;
				return false;
			}
			if (result == 0.0 || math::length(prev - pix) <= eps)
				return true;
			prev = pix;
		}
		return true;
	}
	case line_mode_t::HORIZONTAL:
	{
		if (ipix["y"] - _filter_semi_len - 1 < 0 || ipix["y"] + _filter_semi_len + 1 >= img.height())
			return false;
		img_t row(_filter_semi_len * 2 + 3, 1);
		for (int i = 0; i < _filter_semi_len * 2 + 3; ++i)
			row(i, 0) = img(ipix["x"], ipix["y"] + i - _filter_semi_len - 1);
		math::vec2 temp = { _filter_semi_len + 1, 0 };
		math::vec2 prev = temp;
		double result = 0.0;
		for (int i = 0; i < iteration_count; ++i)
		{
			result = _refine_iteration(row, temp, mode);
			if (result > _max_threshold)
			{
				temp = prev;
				break;
			}
			if (result == 0.0 || math::length(prev - temp) <= eps)
			{
				result = 0.0;
				break;
			}
			prev = temp;
			result = 0.0;
		}
		pix["y"] += temp["x"] - _filter_semi_len - 1;
		return result == 0.0;
	}
	default:
		return false;
	}
}

O3DS::proc_2d::stager_1D_refiner_t::stager_1D_refiner_t() = default;

O3DS::proc_2d::stager_1D_refiner_t::stager_1D_refiner_t(double sigma, double max_threshold) noexcept { init(sigma, max_threshold); }

double O3DS::proc_2d::stager_2D_refiner_t::_refine_iteration(const img_t& img, math::vec2& pix) const noexcept
{
	math::size2 ipix = pix;
	if (ipix["x"] + _filter_semi_len >= img.width() || ipix["x"] < _filter_semi_len
		|| ipix["y"] + _filter_semi_len >= img.height() || ipix["y"] < _filter_semi_len)
		return 0.0;
	math::vec2 T = static_cast<math::vec2>(ipix) - pix;
	filter2D_t gx([this](double x, double y){ return math::d_dt_gauss(x, _sigma) * math::gauss(y, _sigma); },
		{ _filter_semi_len, _filter_semi_len }, { _filter_semi_len, _filter_semi_len }, T);
	filter2D_t gy([this](double x, double y) { return math::gauss(x, _sigma) * math::d_dt_gauss(y, _sigma); },
		{ _filter_semi_len, _filter_semi_len }, { _filter_semi_len, _filter_semi_len }, T);
	filter2D_t gxx([this](double x, double y) { return math::d2_dt2_gauss(x, _sigma) * math::gauss(y, _sigma); },
		{ _filter_semi_len, _filter_semi_len }, { _filter_semi_len, _filter_semi_len }, T);
	filter2D_t gxy([this](double x, double y) { return math::d_dt_gauss(x, _sigma) * math::d_dt_gauss(y, _sigma); },
		{ _filter_semi_len, _filter_semi_len }, { _filter_semi_len, _filter_semi_len }, T);
	filter2D_t gyy([this](double x, double y) { return math::gauss(x, _sigma) * math::d2_dt2_gauss(y, _sigma); },
		{ _filter_semi_len, _filter_semi_len }, { _filter_semi_len, _filter_semi_len }, T);
	double rx = convolution2D(img, gx, ipix);
	double ry = convolution2D(img, gy, ipix);
	double rxx = convolution2D(img, gxx, ipix);
	double rxy = convolution2D(img, gxy, ipix);
	double ryy = convolution2D(img, gyy, ipix);
	math::mat2 H = { rxx, rxy, rxy, ryy };
	math::vec2 l = math::get_eigenvalues(H);
	if (l == math::vec2(math::negative_inf, math::positive_inf))
		return 0.0;
	double lambda = abs(l["x"]) > abs(l["y"]) ? l["x"] : l["y"];
	math::vec2 n = math::get_eigenvector(H, lambda);
	double k = (-(rx * n["x"] + ry * n["y"]) / (rxx * n["x"] * n["x"] + 2.0 * rxy * n["x"] * n["y"] + ryy * n["y"] * n["y"]));
	math::vec2 bias = n * k;
	if (abs(bias["x"]) > 1.0 || abs(bias["y"]) > 1.0)
		return 0.0;
	pix += bias;
	return lambda;
}

void O3DS::proc_2d::stager_2D_refiner_t::init(double sigma, double min_threshold) noexcept
{
	_sigma = sigma;
	_min_threshold = min_threshold;
	_filter_semi_len = 3 * sigma;
}

bool O3DS::proc_2d::stager_2D_refiner_t::refine(const img_t& img, math::vec2& pix, line_mode_t) const noexcept
{
	math::vec2i ipix = pix;
	double result = 0.0;
	math::vec2 prev = pix;
	for (int i = 0; i < iteration_count; ++i)
	{
		result = _refine_iteration(img, pix);
		if (abs(result) < _min_threshold)
		{
			pix = prev;
			return false;
		}
		if (math::length(prev - pix) <= eps)
			return true;
		prev = pix;
	}
	return true;
}

O3DS::proc_2d::stager_2D_refiner_t::stager_2D_refiner_t() = default;

O3DS::proc_2d::stager_2D_refiner_t::stager_2D_refiner_t(double sigma, double min_threshold) noexcept { init(sigma, min_threshold); }

O3DS::proc_2d::rectangle_detector_t::rectangle_detector_t(size_t angle_size, float canny_threshold, size_t line_semi_width, size_t nearest_radius,
	float match_threshold, double equal_threshold) noexcept : _canny_threshold(canny_threshold), _nearest_radius(nearest_radius),
	_match_threshold(match_threshold), _angle_size(angle_size), _equal_threshold(equal_threshold), _line_semi_width(line_semi_width) {}

float O3DS::proc_2d::rectangle_detector_t::get_canny_threshold() const noexcept
{
	return _canny_threshold;
}

void O3DS::proc_2d::rectangle_detector_t::set_canny_threshold(float threshold) noexcept
{
	_canny_threshold = threshold;
}

size_t O3DS::proc_2d::rectangle_detector_t::get_angle_size() const noexcept
{
	return _angle_size;
}

void O3DS::proc_2d::rectangle_detector_t::set_angle_size(size_t angle_size) noexcept
{
	_angle_size = angle_size;
}

float O3DS::proc_2d::rectangle_detector_t::get_match_threshold() const noexcept
{
	return _match_threshold;
}

void O3DS::proc_2d::rectangle_detector_t::set_match_threshold(float threshold) noexcept
{
	_match_threshold = threshold;
}

float O3DS::proc_2d::rectangle_detector_t::get_equal_threshold() const noexcept
{
	return _equal_threshold;
}

void O3DS::proc_2d::rectangle_detector_t::set_equal_threshold(float threshold) noexcept
{
	_equal_threshold = threshold;
}

float O3DS::proc_2d::rectangle_detector_t::get_nearest_radius() const noexcept
{
	return _nearest_radius;
}

void O3DS::proc_2d::rectangle_detector_t::set_nearest_radius(size_t radius) noexcept
{
	_nearest_radius = radius;
}

float O3DS::proc_2d::rectangle_detector_t::get_line_semi_width() const noexcept
{
	return _line_semi_width;
}

void O3DS::proc_2d::rectangle_detector_t::set_line_semi_width(size_t semi_width) noexcept
{
	_line_semi_width = semi_width;
}

void O3DS::proc_2d::rectangle_detector_t::_cpu_find_angles(const O3DS::img_t& temp, uint8_t*& angle_mask) const noexcept
{
	int w = temp.width();
	int h = temp.height();
	angle_mask = new uint8_t[w * h];
	for (int i0 = 0; i0 < h; ++i0)
	{
		for (int j0 = 0; j0 < w; ++j0)
		{
			int counter = 0;
			for (int k = -static_cast<int>(_angle_size); k <= 0; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0 + k, i0) > 0 ? 1 : 0) + (temp(j0, i0 + k) > 0 ? 1 : 0);
			}
			angle_mask[j0 + i0 * w] = std::max(counter > _angle_size ? 1 : 0, (int)angle_mask[j0 + i0 * w]);
			counter = 0;
			for (int k = 0; k <= _angle_size; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0 + k, i0) > 0 ? 1 : 0) + (temp(j0, i0 + k) > 0 ? 1 : 0);
			}
			angle_mask[j0 + i0 * w] = std::max(counter > _angle_size ? 1 : 0, (int)angle_mask[j0 + i0 * w]);
			counter = 0;
			for (int k = -static_cast<int>(_angle_size); k <= 0; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0 + k, i0) > 0 ? 1 : 0);
			}
			for (int k = 0; k <= _angle_size; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0, i0 + k) > 0 ? 1 : 0);
			}
			angle_mask[j0 + i0 * w] = std::max(counter > _angle_size ? 1 : 0, (int)angle_mask[j0 + i0 * w]);
			counter = 0;
			for (int k = 0; k <= _angle_size; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0 + k, i0) > 0 ? 1 : 0);
			}
			for (int k = -static_cast<int>(_angle_size); k <= 0; ++k)
			{
				if (j0 + k <= 0 && j0 + k < w && i0 + k >= 0 && i0 + k < h)
					counter += (temp(j0 + k, i0) > 0 ? 1 : 0);
			}
			angle_mask[j0 + i0 * w] = std::max(counter > _angle_size ? 1 : 0, (int)angle_mask[j0 + i0 * w]);
		}
	}
}

void O3DS::proc_2d::rectangle_detector_t::_gpu_find_angles(const O3DS::img_t& img, uint8_t*& result) const noexcept
{
	static constexpr const char* orientation_kernel_source = R"(
		__kernel void set_orient(__global uchar* in, __global float* out_nx, __global float* out_ny,
			 __global float* filter_xx, __global float* filter_xy, __global float* filter_yy,
			 __global uint* pw, __global uint* ph, __global uint* pfsl)
		{
			int filter_semi_len = pfsl[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				int filter_len = filter_semi_len * 2 + 1;
				float rxx = 0.0;
				for (int i = -filter_semi_len; i <= filter_semi_len; ++i)
				{
					for (int j = -filter_semi_len; j <= filter_semi_len; ++j)
					{
						if (j0 - j < 0 || j0 - j >= w || i0 - i < 0 || i0 - i >= h)
							continue;
						uint zi = j0 - j + (i0 - i) * w;
						uint fy = i + filter_semi_len;
						uint fx = j + filter_semi_len;
						rxx += convert_float(in[zi]) * filter_xx[fx + fy * filter_len];
					}
				}

				float rxy = 0.0;
				for (int i = -filter_semi_len; i <= filter_semi_len; ++i)
				{
					for (int j = -filter_semi_len; j <= filter_semi_len; ++j)
					{
						if (j0 - j < 0 || j0 - j >= w || i0 - i < 0 || i0 - i >= h)
							continue;
						uint zi = j0 - j + (i0 - i) * w;
						uint fy = i + filter_semi_len;
						uint fx = j + filter_semi_len;
						rxy += convert_float(in[zi]) * filter_xy[fx + fy * filter_len];
					}
				}

				float ryy = 0.0;
				for (int i = -filter_semi_len; i <= filter_semi_len; ++i)
				{
					for (int j = -filter_semi_len; j <= filter_semi_len; ++j)
					{
						if (j0 - j < 0 || j0 - j >= w || i0 - i < 0 || i0 - i >= h)
							continue;
						uint zi = j0 - j + (i0 - i) * w;
						uint fy = i + filter_semi_len;
						uint fx = j + filter_semi_len;
						ryy += convert_float(in[zi]) * filter_yy[fx + fy * filter_len];
					}
				}

				float D = (rxx + ryy) * (rxx + ryy) - 4.0f * (rxx * ryy - rxy * rxy);
				if (D >= 0.0f)
				{
					float l1 = (rxx + ryy - sqrt(D)) / 2.0f;
					float l2 = (rxx + ryy + sqrt(D)) / 2.0f;
					float lambda = l1;
					if ((l2 < 0.0f && l1 < 0.0f && lambda > l2) || (l2 > 0.0f && l1 > 0.0f && lambda < l2))
						lambda = l2;
					float nx = 0.0f, ny = 0.0f;
					if (rxx == 0.0f && ryy == 0.0f) {}
					else if (ryy - lambda < 0.000001f && ryy - lambda > -0.000001f)
					{
						ny = 1.0;
						nx = -rxy / (rxx - lambda);
					}
					else
					{
						nx = 1.0;
						ny = -rxy / (ryy - lambda);
					}
					if (nx != 0.0f || ny != 0.0f)
					{
						float l = sqrt(nx * nx + ny * ny);
						nx /= l;
						ny /= l;
						if (lambda < 0.0f)
						{
							nx *= -1.0f;
							ny *= -1.0f;
						}
						out_nx[j0 + i0 * w] = nx;
						out_ny[j0 + i0 * w] = ny;
					}
					else
					{
						out_nx[j0 + i0 * w] = 0.0f;
						out_ny[j0 + i0 * w] = 0.0f;
					}
				}
				else
				{
					out_nx[j0 + i0 * w] = 0.0f;
					out_ny[j0 + i0 * w] = 0.0f;
				}
				if (in[j0 + i0 * w] == 0.0f)
				{
					out_nx[j0 + i0 * w] = 0.0f;
					out_ny[j0 + i0 * w] = 0.0f;
				}
				if (isnan(out_nx[j0 + i0 * w]) || isnan(out_ny[j0 + i0 * w]))
				{
					out_nx[j0 + i0 * w] = 0.0f;
					out_ny[j0 + i0 * w] = 0.0f;
				}
			}
		}
	)";

	static constexpr const char* find_nearest_kernel_source = R"(
		__kernel void linking(__global uchar* in, __global int* out_1, __global int* out_2,
			__global uint* pw, __global uint* ph, __global uint* pr)
		{
			int r = pr[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				int min1 = r * r * 2;
				int bias1 = 0;
				int i1 = 0; int j1 = 0;
				for (int i = -r; i <= r; ++i)
				{
					for (int j = -r; j <= r; ++j)
					{
						int l = j * j + i * i;
						if (l <= min1 && (j != 0 || i != 0) && j0 - j >= 0 && j0 - j < w 
							&& i0 - i >= 0 && i0 - i < h && in[j0 - j + (i0 - i) * w] > 0)
						{
							min1 = l;
							bias1 = (j0 - j + (i0 - i) * w) - (j0 + i0 * w);
							j1 = j;
							i1 = i;
						}
					}
				}
			
				int min2 = r * r * 2;
				int bias2 = 0;
				for (int i = -r; i <= r; ++i)
				{
					for (int j = -r; j <= r; ++j)
					{
						int l = j * j + i * i;
						if (l <= min2 && (j != 0 || i != 0) && j0 - j >= 0 && j0 - j < w 
							&& i0 - i >= 0 && i0 - i < h && in[j0 - j + (i0 - i) * w] > 0
							&& j * j1 + i * i1 <= 0)
						{
							min2 = l;
							bias2 = (j0 - j + (i0 - i) * w) - (j0 + i0 * w);
						}
					}
				}

				if (in[j0 + i0 * w] == 0)
				{
					bias1 = 0;
					bias2 = 0;
				}

				out_1[j0 + i0 * w] = bias1;
				out_2[j0 + i0 * w] = bias2;
			}
		}
	)";

	static constexpr const char* find_angle_kernel_source = R"(
		__kernel void set_angle_mask(__global float* pnx, __global float* pny, __global float* filter,
			__global int* bias_1, __global int* bias_2, __global uchar* out, __global float* rates_out,
			__global uint* pw, __global uint* ph, __global float* pthreshold, __global uint* psize)
		{
			float threshold = pthreshold[0];
			int size = psize[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				float nx0 = pnx[j0 + i0 * w];
				float ny0 = pny[j0 + i0 * w];

				int bias01 = bias_1[j0 + i0 * w];
				int bias02 = bias_2[j0 + i0 * w];

				float nx0b1 = pnx[j0 + i0 * w + bias01];
				float ny0b1 = pny[j0 + i0 * w + bias01];

				float nx0b2 = pnx[j0 + i0 * w + bias02];
				float ny0b2 = pny[j0 + i0 * w + bias02];
	
				if (nx0 > -0.01f && nx0 < 0.01f && ny0 > -0.01f && ny0 < 0.01f
					&& !(nx0b1 > -0.01f && nx0b1 < 0.01f && ny0b1 > -0.01f && ny0b1 < 0.01f)
					&& !(nx0b2 > -0.01f && nx0b2 < 0.01f && ny0b2 > -0.01f && ny0b2 < 0.01f))
				{
					nx0 = (nx0b1 + nx0b2) / 2.0f;
					ny0 = (ny0b1 + ny0b2) / 2.0f;
				}

				if (nx0 > -0.01f && nx0 < 0.01f && ny0 > -0.01f && ny0 < 0.01f)
				{
					out[j0 + i0 * w] = 0;
					rates_out[j0 + i0 * w] = 0.0f;
				}
				else
				{
					float prev_nx = nx0;
					float prev_ny = ny0;
					int prev_index = j0 + i0 * w;
					int b_complete = 1;
					for (int i = 0; i < size; ++i)
					{
						int prev_i = prev_index / w;
						int prev_j = prev_index - prev_i * w;
						int bias1 = bias_1[prev_j + prev_i * w];
						int bias2 = bias_2[prev_j + prev_i * w];
						if (bias1 != 0 && bias2 != 0)
						{
							float nx1 = pnx[prev_j + prev_i * w + bias1];
							float ny1 = pny[prev_j + prev_i * w + bias1];
							if (nx1 * prev_nx + ny1 * prev_ny < 0.0f)
							{
								nx1 = -nx1;
								ny1 = -ny1;
							}
							float nx2 = pnx[prev_j + prev_i * w + bias2];
							float ny2 = pny[prev_j + prev_i * w + bias2];
							if (nx2 * prev_nx + ny2 * prev_ny < 0.0f)
							{
								nx2 = -nx2;
								ny2 = -ny2;
							}
							int index1 = prev_index + bias1;
							int cur_i1 = (index1 / w) - prev_i;
							int cur_j1 = (index1 - (index1 / w) * w) - prev_j;
							int index2 = prev_index + bias2;
							int cur_i2 = (index2 / w) - prev_i;
							int cur_j2 = (index2 - (index2 / w) * w) - prev_j;
							float px = -prev_ny;
							float py = prev_nx;
							if (px * cur_j1 + py * cur_i1 > 0)
							{
								prev_index = index1;
								prev_nx = nx1;
								prev_ny = ny1;
							}
							else
							{
								prev_index = index2;
								prev_nx = nx2;
								prev_ny = ny2;
							}
						}
						else
							b_complete = 0;
					}
					float integral = 0.0f;
					int beg_index = prev_index;
					if (b_complete)
					{
						for (int i = -size; i <= size; ++i)
						{
							int prev_i = prev_index / w;
							int prev_j = prev_index - prev_i * w;
							int bias1 = bias_1[prev_j + prev_i * w];
							int bias2 = bias_2[prev_j + prev_i * w];
							if (bias1 != 0 && bias2 != 0)
							{
								float nx1 = pnx[prev_j + prev_i * w + bias1];
								float ny1 = pny[prev_j + prev_i * w + bias1];
								if (nx1 * prev_nx + ny1 * prev_ny < 0.0f)
								{
									nx1 = -nx1;
									ny1 = -ny1;
								}
								float nx2 = pnx[prev_j + prev_i * w + bias2];
								float ny2 = pny[prev_j + prev_i * w + bias2];
								if (nx2 * prev_nx + ny2 * prev_ny < 0.0f)
								{
									nx2 = -nx2;
									ny2 = -ny2;
								}
								int index1 = prev_index + bias1;
								int cur_i1 = (index1 / w) - prev_i;
								int cur_j1 = (index1 - (index1 / w) * w) - prev_j;
								int index2 = prev_index + bias2;
								int cur_i2 = (index2 / w) - prev_i;
								int cur_j2 = (index2 - (index2 / w) * w) - prev_j;
								float px = prev_ny;
								float py = -prev_nx;
								if (px * cur_j1 + py * cur_i1 > 0)
								{
									prev_index = index1;
									float ndot = 1.0f - prev_nx * nx1 - prev_ny * ny1;
									ndot = ndot < 0.0f ? -ndot : ndot;
									ndot = ndot < 0.001 ? 0.001 : ndot;
									prev_nx = nx1;
									prev_ny = ny1;
									float f = filter[i + size];
									integral += f * ndot;
								}
								else
								{
									prev_index = index2;
									float ndot = 1.0f - prev_nx * nx2 - prev_ny * ny2;
									ndot = ndot < 0.0f ? -ndot : ndot;
									ndot = ndot < 0.001 ? 0.001 : ndot;
									prev_nx = nx2;
									prev_ny = ny2;
									float f = filter[i + size];
									integral += f * ndot;
								}
							}
						}
						int i1 = beg_index / w;
						int j1 = beg_index - i1 * w;
						int i2 = prev_index / w;
						int j2 = prev_index - i2 * w;
						float px1 = j1 - j0;
						float py1 = i1 - i0;
						float l1 = sqrt(px1 * px1 + py1 * py1);
						px1 /= l1;
						py1 /= l1;
						float px2 = j2 - j0;
						float py2 = i2 - i0;
						float l2 = sqrt(px2 * px2 + py2 * py2);
						px2 /= l2;
						py2 /= l2;
						float res_dot = px1 * px2 + py1 * py2;
						res_dot = res_dot < 0.0f ? -res_dot : res_dot;
						if (res_dot > 0.7f || l1 < (size / 2) || l2 < (size / 2))
							integral = 0.0f;
						if (integral > threshold)
							out[j0 + i0 * w] = 1;
						else
							out[j0 + i0 * w] = 0;
						rates_out[j0 + i0 * w] = integral;
					}
					else
						out[j0 + i0 * w] = 0;
				}
			}
		}
	)";

	static constexpr const char* twins_detection_kernel_source = R"(
		__kernel void no_twins(__global uchar* mask, __global float* rates,
			__global uchar* out_mask, __global uint* pw, __global uint* ph, __global uint* psize)
		{
			int size = psize[0];
			int w = pw[0];
			int h = ph[0];

			int index = get_global_id(0);
			int i0 = convert_float(index) / convert_float(w);
			int j0 = index - (i0 * w);

			if (index < w * h)
			{
				out_mask[j0 + i0 * w] = mask[j0 + i0 * w];			

				if (mask[j0 + i0 * w] > 0)
				{
					for (int i = -size; i <= size; ++i)
					{
						for (int j = -size; j <= size; ++j)
						{
							if (mask[j0 - j + (i0 - i) * w] > 0 
								&& rates[j0 - j + (i0 - i) * w] > rates[j0 + i0 * w])
								out_mask[j0 + i0 * w] = 0;
						}
					}
				}
			}
		}
	)";

	boost::compute::context gpu_context(gpu_device);
	boost::compute::command_queue gpu_command_queue(gpu_context, gpu_device);

	constexpr cl_uint filter_semi_len = 10;
	constexpr cl_uint filter_len = filter_semi_len * 2 + 1;
	constexpr float sigma = 1.5;

	filter1D_t classic_gauss([&sigma](double x) { return 100.0 * math::gauss(x, sigma) - 0.25; }, _angle_size, _angle_size, 0);
	filter2D_t gxx([sigma, filter_semi_len](double x, double y) { return math::d2_dt2_gauss(x, sigma) * math::gauss(y, sigma); },
		{ filter_semi_len, filter_semi_len }, { filter_semi_len, filter_semi_len }, 0.0f);
	filter2D_t gxy([sigma, filter_semi_len](double x, double y) { return math::d_dt_gauss(x, sigma) * math::d_dt_gauss(y, sigma); },
		{ filter_semi_len, filter_semi_len }, { filter_semi_len, filter_semi_len }, 0.0f);
	filter2D_t gyy([sigma, filter_semi_len](double x, double y) { return math::gauss(x, sigma) * math::d2_dt2_gauss(y, sigma); },
		{ filter_semi_len, filter_semi_len }, { filter_semi_len, filter_semi_len }, 0.0f);

	// Define orientation of curves
	const uint8_t* input_data = img.get_format()->bitmap;
	const float* filter_xx_data = gxx.data();
	const float* filter_xy_data = gxy.data();
	const float* filter_yy_data = gyy.data();
	cl_uint w = img.width();
	cl_uint h = img.height();
	boost::compute::program program1 = boost::compute::program::create_with_source(orientation_kernel_source, gpu_context);
	program1.build();
	boost::compute::kernel kernel1(program1, "set_orient");
	boost::compute::buffer input(gpu_context, w * h);
	boost::compute::buffer out_nx(gpu_context, w * h * sizeof(cl_float));
	boost::compute::buffer out_ny(gpu_context, w * h * sizeof(cl_float));
	boost::compute::buffer filter_xx(gpu_context, filter_len * filter_len * sizeof(cl_float));
	boost::compute::buffer filter_xy(gpu_context, filter_len * filter_len * sizeof(cl_float));
	boost::compute::buffer filter_yy(gpu_context, filter_len * filter_len * sizeof(cl_float));
	boost::compute::buffer wb(gpu_context, sizeof(cl_uint));
	boost::compute::buffer hb(gpu_context, sizeof(cl_uint));
	boost::compute::buffer fslb(gpu_context, sizeof(cl_uint));
	gpu_command_queue.enqueue_write_buffer(input, 0, w * h, input_data);
	gpu_command_queue.enqueue_write_buffer(filter_xx, 0, filter_len* filter_len * sizeof(cl_float), filter_xx_data);
	gpu_command_queue.enqueue_write_buffer(filter_xy, 0, filter_len* filter_len * sizeof(cl_float), filter_xy_data);
	gpu_command_queue.enqueue_write_buffer(filter_yy, 0, filter_len* filter_len * sizeof(cl_float), filter_yy_data);
	gpu_command_queue.enqueue_write_buffer(wb, 0, sizeof(cl_uint), &w);
	gpu_command_queue.enqueue_write_buffer(hb, 0, sizeof(cl_uint), &h);
	gpu_command_queue.enqueue_write_buffer(fslb, 0, sizeof(cl_uint), &filter_semi_len);
	kernel1.set_arg(0, input);
	kernel1.set_arg(1, out_nx);
	kernel1.set_arg(2, out_ny);
	kernel1.set_arg(3, filter_xx);
	kernel1.set_arg(4, filter_xy);
	kernel1.set_arg(5, filter_yy);
	kernel1.set_arg(6, wb);
	kernel1.set_arg(7, hb);
	kernel1.set_arg(8, fslb);
	gpu_command_queue.enqueue_1d_range_kernel(kernel1, 0, get_range_size(w * h, 128), 128);

	// Define curves as linked points
	boost::compute::program program2 = boost::compute::program::create_with_source(find_nearest_kernel_source, gpu_context);
	program2.build();
	boost::compute::kernel kernel2(program2, "linking");
	unsigned nearest_radius = _nearest_radius;
	boost::compute::buffer nearest_1(gpu_context, w * h * sizeof(cl_int));
	boost::compute::buffer nearest_2(gpu_context, w * h * sizeof(cl_int));
	boost::compute::buffer rb(gpu_context, sizeof(cl_uint));
	gpu_command_queue.enqueue_write_buffer(rb, 0, sizeof(cl_uint), &nearest_radius);
	kernel2.set_arg(0, input);
	kernel2.set_arg(1, nearest_1);
	kernel2.set_arg(2, nearest_2);
	kernel2.set_arg(3, wb);
	kernel2.set_arg(4, hb);
	kernel2.set_arg(5, rb);
	gpu_command_queue.enqueue_1d_range_kernel(kernel2, 0, get_range_size(w * h, 128), 128);

	// Define angles candidates
	boost::compute::program program3 = boost::compute::program::create_with_source(find_angle_kernel_source, gpu_context);
	program3.build();
	boost::compute::kernel kernel3(program3, "set_angle_mask");
	unsigned angle_size = _angle_size;
	const float* classic_filter_data = classic_gauss.data();
	boost::compute::buffer rates_out(gpu_context, w * h * sizeof(float));
	boost::compute::buffer angle_filter(gpu_context, (2 * _angle_size + 1) * (2 * _angle_size + 1) * sizeof(cl_float));
	boost::compute::buffer thrb(gpu_context, sizeof(cl_float));
	boost::compute::buffer asizeb(gpu_context, sizeof(cl_uint));
	gpu_command_queue.enqueue_write_buffer(thrb, 0, sizeof(cl_float), &_match_threshold);
	gpu_command_queue.enqueue_write_buffer(asizeb, 0, sizeof(cl_float), &angle_size);
	gpu_command_queue.enqueue_write_buffer(angle_filter, 0, (2 * _angle_size + 1)* (2 * _angle_size + 1) * sizeof(cl_float), classic_filter_data);
	kernel3.set_arg(0, out_nx);
	kernel3.set_arg(1, out_ny);
	kernel3.set_arg(2, angle_filter);
	kernel3.set_arg(3, nearest_1);
	kernel3.set_arg(4, nearest_2);
	kernel3.set_arg(5, input);
	kernel3.set_arg(6, rates_out);
	kernel3.set_arg(7, wb);
	kernel3.set_arg(8, hb);
	kernel3.set_arg(9, thrb);
	kernel3.set_arg(10, asizeb);
	gpu_command_queue.enqueue_1d_range_kernel(kernel3, 0, get_range_size(w * h, 128), 128);

	// Angles filtration
	result = new uint8_t[w * h];
	boost::compute::program program4 = boost::compute::program::create_with_source(twins_detection_kernel_source, gpu_context);
	program4.build();
	boost::compute::kernel kernel4(program4, "no_twins");
	boost::compute::buffer output(gpu_context, w * h);
	kernel4.set_arg(0, input);
	kernel4.set_arg(1, rates_out);
	kernel4.set_arg(2, output);
	kernel4.set_arg(3, wb);
	kernel4.set_arg(4, hb);
	kernel4.set_arg(5, asizeb);
	gpu_command_queue.enqueue_1d_range_kernel(kernel4, 0, get_range_size(w * h, 128), 128);
	gpu_command_queue.enqueue_read_buffer(input, 0, w * h, result);
}

void O3DS::proc_2d::rectangle_detector_t::get_rectangles(const img_t& img, std::vector<rect_t>& rects) const noexcept
{
	img_t temp = canny_edge_detection(img, 1.7, _canny_threshold);
	int w = temp.width();
	int h = temp.height();
	uint8_t* angle_mask;
	if (compute::use_gpu)
		_gpu_find_angles(temp, angle_mask);
	else
		_cpu_find_angles(temp, angle_mask);

	std::vector<math::vec2> angles;
	angles.reserve(w * h / 8);
	stager_2D_refiner_t refiner(1.5);
	line_detector_t line_detector(&refiner);

#pragma omp parallel for
	for (int index = 0; index < w * h; ++index)
	{
		int i = index / w;
		int j = index - i * w;
		if (angle_mask[j + i * w])
		{
			math::vec2 angle = line_detector.get_point(temp, { j, i });
#pragma omp critical
			{
				angles.push_back(angle);
			}
		}
	}

	std::vector<char> true_mask(angles.size());
	for (auto i = angles.begin(); i != angles.end(); ++i)
	{
		true_mask[std::distance(angles.begin(), i)] = 1;
		for (auto j = i + 1; j != angles.end(); ++j)
		{
			if (math::length(*j - *i) <= _angle_size * 2 + 1)
			{
				true_mask[std::distance(angles.begin(), i)] = 0;
				break;
			}
		}
	}

	delete[] angle_mask;
	int angle_count = angles.size();
	std::vector<std::tuple<math::vec2, math::vec2>> lines;
	lines.reserve(angle_count * angle_count);

#pragma omp parallel for
	for (int index = 0; index < angle_count * angle_count; ++index)
	{
		int i = index / angle_count;
		int j = index - i * angle_count;
		if (j > i && true_mask[i] && true_mask[j])
		{
			math::vec2 a = angles[i];
			math::vec2 b = angles[j];
			int h_count = 20;
			math::vec2 p = (b - a);
			if (math::length(p) < h_count)
				h_count = 0;
			p = p * (1.0 / static_cast<double>(h_count));
			math::vec2 n = { -p["y"], p["x"] };
			h_count = math::length(p) < 0.5 ? 0 : h_count;
			math::vec2 current = a;
			int white_counter = 0; int black_counter = 0;
			for (int k = 0; k < h_count; ++k)
			{
				current += p;
				math::size2 pix = current;
				int radius = _line_semi_width;
				if (pix["x"] > radius && pix["y"] > radius 
					&& pix["x"] < w - (radius + 1) && pix["y"] < h - (radius + 1))
				{
					++black_counter;
					for (int t = -radius; t <= radius; ++t)
					{
						math::size2 x;
						if (abs(n["x"]) > abs(n["y"])) x = { pix["x"] + t, pix["y"] };
						else x = { pix["x"], pix["y"] + t };
						if (temp(x["x"], x["y"]) > 0)
						{
							--black_counter;
							++white_counter;
							break;
						}
					}
				}
			}
			if (black_counter <= 1)
			{
#pragma omp critical
				{
					lines.emplace_back(a, b);
				}
			}
		}
	}

	int line_count = lines.size();
	std::vector<std::pair<math::vec2, math::vec2>> true_lines;
	true_lines.reserve(line_count);

#pragma omp parallel for
	for (int i = 0; i < line_count; ++i)
	{
		math::vec2 p1 = std::get<1>(lines[i]) - std::get<0>(lines[i]);
		math::vec2 a1 = std::get<0>(lines[i]);
		math::vec2 b1 = std::get<1>(lines[i]);
		math::vec2 n = { -p1["y"], p1["x"] };
		bool b_emplace = true;
		for (int j = 0; j < line_count; ++j)
		{
			if (i == j) continue;
			math::vec2 p2 = std::get<1>(lines[j]) - std::get<0>(lines[j]);
			math::vec2 a2 = std::get<0>(lines[j]);
			math::vec2 b2 = std::get<1>(lines[j]);
			double lna = abs(math::dot(math::normalize(p1), a1 - a2));
			lna = sqrt(math::length(a1 - a2) * math::length(a1 - a2) - lna * lna);
			lna = isnan(lna) ? 0.0 : lna;
			double lnb = abs(math::dot(math::normalize(p1), b1 - b2));
			lna = sqrt(math::length(b1 - b2) * math::length(b1 - b2) - lnb * lnb);
			lnb = isnan(lnb) ? 0.0 : lnb;
			double ln = (lna + lnb) / 2.0;
			double lp1 = math::length(p1);
			double lp2 = math::length(p2);
			math::vec2 avg = (a1 + b1 + a2 + b2) * 0.25;
			double ndot1 = math::dot(a1 - avg, b1 - avg);
			double ndot2 = math::dot(a2 - avg, b2 - avg);
			bool bdot = (ndot1 < 0.0 || ndot2 < 0.0);
			if (abs(math::dot(math::normalize(p1), math::normalize(p2))) > 0.80 && ln < 2.0 && bdot && lp2 > lp1)
			{
				b_emplace = false;
				break;
			}
		}
		if (b_emplace)
		{
#pragma omp critical
			{
				true_lines.emplace_back(std::get<0>(lines[i]), std::get<1>(lines[i]));
			}
		}
	}

	line_count = true_lines.size();
	std::vector<std::tuple<math::vec2, math::vec2, math::vec2>> triples;
	triples.reserve(line_count * line_count);
	for (int index = 0; index < line_count * line_count; ++index)
	{
		int i = index / line_count;
		int j = index - i * line_count;
		if (j > i && true_lines[i].first == true_lines[j].first)
			triples.emplace_back(true_lines[i].second, true_lines[i].first, true_lines[j].second);
		else if (j > i && true_lines[i].first == true_lines[j].second)
			triples.emplace_back(true_lines[i].second, true_lines[i].first, true_lines[j].first);
		else if (j > i && true_lines[i].second == true_lines[j].first)
			triples.emplace_back(true_lines[i].first, true_lines[i].second, true_lines[j].second);
		else if (j > i && true_lines[i].second == true_lines[j].second)
			triples.emplace_back(true_lines[i].first, true_lines[i].second, true_lines[j].first);
	}

	int triple_count = triples.size();
	std::vector<std::pair<bool, rect_t>> temp_rects;
	temp_rects.reserve(triple_count);
	for (int i = 0; i < triple_count; ++i)
	{
		rect_t rect;
		math::vec2 a1 = std::get<0>(triples[i]);
		math::vec2 a2 = std::get<1>(triples[i]);
		math::vec2 a3 = std::get<2>(triples[i]);
		math::vec2 a4 = { -1.0, -1.0 };
		for (int j = i + 1; j < triple_count; ++j)
		{
			math::vec2 b1 = std::get<0>(triples[j]);
			math::vec2 b2 = std::get<1>(triples[j]);
			math::vec2 b3 = std::get<2>(triples[j]);
			if ((b1 == a1 && b3 == a3) || (b1 == a3 && b3 == a1))
			{ a4 = b2; break; }
		}
		if (a4["x"] < 0.0 || a4["y"] < 0.0)
		{
			for (int j = i + 1; j < triple_count; ++j)
			{
				math::vec2 b1 = std::get<0>(triples[j]);
				math::vec2 b2 = std::get<1>(triples[j]);
				math::vec2 b3 = std::get<2>(triples[j]);
				if (b1 == a2 && b2 == a3) { a4 = b3; break; }
				if (b3 == a2 && b2 == a3) { a4 = b1; break; }
				if (b1 == a2 && b2 == a1) { a4 = b3; break; }
				if (b3 == a2 && b2 == a1) { a4 = b1; break; }
			}
		}
		if (a4["x"] >= 0.0 && a4["y"] >= 0.0)
		{
			rect.verts[0] = a1; rect.verts[1] = a2;
			rect.verts[2] = a3; rect.verts[3] = a4;
			temp_rects.emplace_back(true, rect);
		}
	}

	rects.reserve(temp_rects.size());
	for (int i = 0; i < temp_rects.size(); ++i)
	{
		for (int j = i + 1; j < temp_rects.size() && temp_rects[i].first; ++j)
		{
			if (rect_t::equal(temp_rects[i].second, temp_rects[j].second, _equal_threshold))
				temp_rects[j].first = false;
		}
		if (temp_rects[i].first) rects.push_back(temp_rects[i].second);
	}
}

double O3DS::proc_2d::rectangle_detector_t::rect_t::square() const noexcept
{
	math::vec2 diag1 = verts[2] - verts[0];
	math::vec2 diag2 = verts[3] - verts[1];
	return math::length(diag1) * math::length(diag2) / 2.0;
}

bool O3DS::proc_2d::rectangle_detector_t::rect_t::equal(const rect_t& a, const rect_t& b, double threshold) noexcept
{
	math::vec2 diag1a = a.verts[0] - a.verts[2];
	math::vec2 diag2a = a.verts[1] - a.verts[3];
	math::vec2 diag1b = b.verts[0] - b.verts[2];
	math::vec2 diag2b = b.verts[1] - b.verts[3];
	double dota = abs(math::dot(diag1a, diag2a));
	double dotb = abs(math::dot(diag1b, diag2b));
	return abs(dota - dotb) < threshold;
}
