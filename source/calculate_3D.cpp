#include "calculate_3D.hpp"
#include "image.hpp"
#include "numerics.hpp"
#include <vector>
#include <list>
#include <algorithm>

void O3DS::proc_3d::cameras_pair_t::_calc_sides(const std::array<math::vec2, 4>& points,
	double& out_a, double& out_b) noexcept
{
	math::vec2 center = (points[0] + points[1] + points[2] + points[3]) * 0.25;
	math::vec2 i = { 1.0, 0.0 }, j = { 0.0, 1.0 };
	std::pair<math::vec2, math::vec2> horizontal[2], vertical[2];
	for (const auto& a : points)
	{
		for (const auto& b : points)
		{
			if (a == b) continue;
			if (math::dot(a - center, j) > 0 && math::dot(b - center, j) > 0)
				horizontal[0].first = a, horizontal[0].second = b;
			else if (math::dot(a - center, i) > 0 && math::dot(b - center, i) > 0)
				vertical[0].first = a, vertical[0].second = b;
			else if (math::dot(a - center, j) < 0 && math::dot(b - center, j) < 0)
				horizontal[1].first = a, horizontal[1].second = b;
			else if (math::dot(a - center, i) < 0 && math::dot(b - center, i) < 0)
				vertical[1].first = a, vertical[1].second = b;
		}
	}

	double a0 = math::length(horizontal[0].first - horizontal[0].second);
	double a1 = math::length(horizontal[1].first - horizontal[1].second);
	out_a = (a0 + a1) / 2.0;
	double b0 = math::length(vertical[0].first - vertical[0].second);
	double b1 = math::length(vertical[1].first - vertical[1].second);
	out_b = (b0 + b1) / 2.0;
}

bool O3DS::proc_3d::cameras_pair_t::_calc_calibration_params(const img_t& source, double width, double height, double depth,
	const proc_2d::rectangle_detector_t& detector, const math::vec2 pixel_size, double& out_h) noexcept
{
	std::vector<proc_2d::rectangle_detector_t::rect_t> rects;
	detector.get_rectangles(source, rects);
	if (rects.size() != 2) return false;
	struct { double a, b; } rect_0, rect_1;
	_calc_sides(rects[0].verts, rect_0.a, rect_0.b);
	_calc_sides(rects[1].verts, rect_1.a, rect_1.b);
	if (rect_0.a * rect_0.b > rect_1.a * rect_1.b)
	{
		decltype(rect_0) temp_rect = rect_0;
		rect_0 = rect_1;
		rect_1 = temp_rect;
	}
	double h_hor = depth / (width / (2.0 * rect_0.a * pixel_size["x"]) - width / (2.0 * rect_1.a * pixel_size["x"]));
	double h_vert = depth / (height / (2.0 * rect_0.b * pixel_size["y"]) - height / (2.0 * rect_1.b * pixel_size["y"]));
	out_h = (h_hor + h_vert) / 2.0;
	return true;
}

O3DS::proc_3d::cameras_pair_t::cameras_pair_t(const camera_data_t& cam_a, const camera_data_t& cam_b,
	const calibration_data_t& data) noexcept
{
	if (_calc_calibration_params(cam_a.source, data.width, data.height, data.depth,
			cam_a.detector, cam_a.pixel_size, _cama_h)
		&& _calc_calibration_params(cam_b.source, data.width, data.height, data.depth,
			cam_b.detector, cam_b.pixel_size, _camb_h))
		_calibration_status = true;
	else
		_calibration_status = false;
	_cama_size = { cam_a.source.width(), cam_a.source.height() };
	_camb_size = { cam_b.source.width(), cam_b.source.height() };
	_cama_pixel_size = cam_a.pixel_size;
	_camb_pixel_size = cam_b.pixel_size;
}

O3DS::proc_3d::cameras_pair_t::cameras_pair_t(double cama_h, math::size2 cama_size, double camb_h, math::size2 camb_size) noexcept
	: _cama_h(cama_h), _cama_size(cama_size), _camb_h(camb_h), _camb_size(camb_size), _calibration_status(true) {}

bool O3DS::proc_3d::cameras_pair_t::get_calibration_status() const noexcept
{
	return _calibration_status;
}

double O3DS::proc_3d::cameras_pair_t::get_camera_h(camera_selection_t pick) const noexcept
{
	return pick == camera_selection_t::A ? _cama_h : _camb_h;
}

void O3DS::proc_3d::cameras_pair_t::get_camera_rotation(proc_2d::rectangle_detector_t::rect_t source,
	math::vec3& target, camera_selection_t pick) const noexcept
{
	std::array<math::vec2, 4>& verts = source.verts;
	math::vec2 a1s = verts[0]; math::vec2 b1s = verts[2];
	math::vec2 ps = b1s - a1s;
	math::vec2 a2s = verts[1]; math::vec2 b2s = verts[3];
	math::vec2 qs = b2s - a2s;
	double ts = (ps["x"] * (a2s["y"] - a1s["y"]) + ps["y"] * (a1s["x"] - a2s["x"])) / (ps["y"] * qs["x"] - ps["x"] * qs["y"]);
	math::vec2 source_center = a2s + qs * ts;
	math::vec2 img_center = pick == camera_selection_t::A ? static_cast<math::vec2>(_cama_size) * 0.5
		: static_cast<math::vec2>(_camb_size) * 0.5;
	for (auto& vert : verts)
	{
		vert -= img_center;
		vert["x"] *= pick == camera_selection_t::A ? _cama_pixel_size["x"] : _camb_pixel_size["x"];
		vert["y"] *= pick == camera_selection_t::A ? _cama_pixel_size["y"] : _camb_pixel_size["y"];
	}
	source_center -= img_center;
	source_center["x"] *= pick == camera_selection_t::A ? _cama_pixel_size["x"] : _camb_pixel_size["x"];
	source_center["y"] *= pick == camera_selection_t::A ? _cama_pixel_size["y"] : _camb_pixel_size["y"];

	int left_top_num = std::distance(verts.begin(), std::max_element(verts.begin(), verts.end(),
		[](const math::vec2& a, const math::vec2& b) { return math::dot(a, math::vec2(-1.0, -1.0)) < math::dot(b, math::vec2(-1.0, -1.0)); }));
	int left_buttom_num = std::distance(verts.begin(), std::max_element(verts.begin(), verts.end(),
		[](const math::vec2& a, const math::vec2& b) { return math::dot(a, math::vec2(-1.0, 1.0)) < math::dot(b, math::vec2(-1.0, 1.0)); }));
	int right_top_num = std::distance(verts.begin(), std::max_element(verts.begin(), verts.end(),
		[](const math::vec2& a, const math::vec2& b) { return math::dot(a, math::vec2(1.0, -1.0)) < math::dot(b, math::vec2(1.0, -1.0)); }));
	int right_buttom_num = std::distance(verts.begin(), std::max_element(verts.begin(), verts.end(),
		[](const math::vec2& a, const math::vec2& b) { return math::dot(a, math::vec2(1.0, 1.0)) < math::dot(b, math::vec2(1.0, 1.0)); }));

	std::array<double, 4> lens; std::array<int, 4> order;
	lens[0] = math::length(verts[left_top_num] - source_center);
	lens[2] = math::length(verts[right_buttom_num] - source_center);
	lens[1] = math::length(verts[right_top_num] - source_center);
	lens[3] = math::length(verts[left_buttom_num] - source_center);
	int max_num = std::distance(lens.begin(), std::max_element(lens.begin(), lens.end()));
	order[0] = max_num; lens[max_num] = 0;
	max_num = max_num + 2 <= 3 ? max_num + 2 : max_num - 2;
	order[2] = max_num; lens[max_num] = 0;
	max_num = std::distance(lens.begin(), std::max_element(lens.begin(), lens.end()));
	order[1] = max_num; lens[max_num] = 0;
	max_num = max_num + 2 <= 3 ? max_num + 2 : max_num - 2;
	order[3] = max_num; lens[max_num] = 0;

	double nc = source_center["x"], mc = source_center["y"];
	double xc = 0.0, yc = 0.0, zc = 0.0;
	double h = pick == camera_selection_t::A ? _cama_h : _camb_h;
	double n1 = verts[order[0]]["x"], n2 = verts[order[1]]["x"], n3 = verts[order[2]]["x"], n4 = verts[order[3]]["x"];
	double m1 = verts[order[0]]["y"], m2 = verts[order[1]]["y"], m3 = verts[order[2]]["y"], m4 = verts[order[3]]["y"];
	math::matrix<double, 6, 6> koeffs = {
		1.0, 0.0, 0.0, 0.0, -n3 / 2.0, 0.0,
		0.0, 1.0, 0.0, 0.0, -m3 / 2.0, 0.0,
		0.0, 0.0, 1.0, 0.0, -h / 2.0, 0.0,
		1.0, 0.0, 0.0, -n2 / 2.0, 0.0, -n4 / 2.0,
		0.0, 1.0, 0.0, -m2 / 2.0, 0.0, -m4 / 2.0,
		0.0, 0.0, 1.0, -h / 2.0, 0.0, -h / 2.0
	};
	math::vector<double, 6> free_terms = { n1 / 2.0, m1 / 2.0, -h / 2.0, 0.0, 0.0, -h };
	math::vector<double, 6> solution = math::get_slau_solution(koeffs, free_terms);
	xc = solution[0]; yc = solution[1]; zc = solution[2];
	double t1 = 1.0, t2 = solution[3], t3 = solution[4], t4 = solution[5];

	std::array<double, 4> t = { t1, t2, t3, t4 }, n = { n1, n2, n3, n4 }, m = { m1, m2, m3, m4 };
	std::array<math::vec3, 4> target_verts;
	int target_left_top_num = 0, target_left_buttom_num = 0, target_right_top_num = 0, target_right_buttom_num = 0;
	for (int i = 0; i < 4; ++i)
	{
		if (order[i] == left_top_num)
			target_left_top_num = i;
		else if (order[i] == left_buttom_num)
			target_left_buttom_num = i;
		else if (order[i] == right_top_num)
			target_right_top_num = i;
		else if (order[i] == right_buttom_num)
			target_right_buttom_num = i;
	}
	target_verts[target_left_top_num] = math::vec3(0.0, 0.0, -h)
		+ (math::vec3(n[target_left_top_num], m[target_left_top_num], 0.0) - math::vec3(0.0, 0.0, -h)) * t[target_left_top_num];
	target_verts[target_left_buttom_num] = math::vec3(0.0, 0.0, -h)
		+ (math::vec3(n[target_left_buttom_num], m[target_left_buttom_num], 0.0) - math::vec3(0.0, 0.0, -h)) * t[target_left_buttom_num];
	target_verts[target_right_top_num] = math::vec3(0.0, 0.0, -h)
		+ (math::vec3(n[target_right_top_num], m[target_right_top_num], 0.0) - math::vec3(0.0, 0.0, -h)) * t[target_right_top_num];
	target_verts[target_right_buttom_num] = math::vec3(0.0, 0.0, -h)
		+ (math::vec3(n[target_right_buttom_num], m[target_right_buttom_num], 0.0) - math::vec3(0.0, 0.0, -h)) * t[target_right_buttom_num];

	math::vec3 i = math::normalize(target_verts[target_right_top_num] - target_verts[target_left_top_num]);
	math::vec3 j = math::normalize(target_verts[target_left_buttom_num] - target_verts[target_left_top_num]);
	math::vec3 k = math::normalize(math::cross(i, j));
	math::vec3 i0 = math::normalize(math::vec3(i["x"], 0.0, i["y"]));
	math::vec3 k0 = math::normalize(math::vec3(k["x"], 0.0, k["z"]));
	double psi = acos(math::dot(k0, math::vec3(0.0, 0.0, 1.0)));
	psi = k0["x"] > 0.0 ? -psi : psi;
	double tetta = acos(math::dot(k0, k));
	tetta = k0["y"] > 0.0 ? -tetta : tetta;
	double gamma = acos(math::dot(i0, i));
	gamma = i0["y"] > 0.0 ? -gamma : gamma;
	target = { psi, tetta, gamma };
}

void O3DS::proc_3d::cameras_pair_t::get_cameras_position(std::vector<proc_2d::rectangle_detector_t::rect_t> rects_a,
	std::vector<proc_2d::rectangle_detector_t::rect_t> rects_b,
	math::vec3 cama_euler, math::vec3 camb_euler, math::vec3& cama_position, math::vec3& camb_position) noexcept
{
	if (rects_a.size() != rects_b.size() || rects_a.size() < 6)
		return;

	std::vector<math::vec2> markers_a(rects_a.size()), markers_b(rects_b.size());
	math::vec2 cama_center = static_cast<math::vec2>(_cama_size) * 0.5;
	for (int i = 0; i < rects_a.size(); ++i)
	{
		for (auto& vert : rects_a[i].verts)
		{
			vert -= cama_center;
			vert["x"] *= _cama_pixel_size["x"];
			vert["y"] *= _cama_pixel_size["y"];
		}
		markers_a[i] = (rects_a[i].verts[0] + rects_a[i].verts[1] + rects_a[i].verts[2] + rects_a[i].verts[3]) * 0.25;
	}
	math::vec2 camb_center = static_cast<math::vec2>(_camb_size) * 0.5;
	for (int i = 0; i < rects_b.size(); ++i)
	{
		for (auto& vert : rects_b[i].verts)
		{
			vert -= camb_center;
			vert["x"] *= _camb_pixel_size["x"];
			vert["y"] *= _camb_pixel_size["y"];
		}
		markers_b[i] = (rects_b[i].verts[0] + rects_b[i].verts[1] + rects_b[i].verts[2] + rects_b[i].verts[3]) * 0.25;
	}

	math::mat3 dcm_a = math::dcm_12g(cama_euler);
	math::mat3 dcm_b = math::dcm_12g(camb_euler);
	std::vector<double> a1 = { markers_a[0]["x"], markers_a[1]["x"], markers_a[2]["x"] };
	std::vector<double> b1 = { markers_a[0]["y"], markers_a[1]["y"], markers_a[2]["y"] };
	std::vector<double> a2 = { markers_b[0]["x"], markers_b[1]["x"], markers_b[2]["x"] };
	std::vector<double> b2 = { markers_b[0]["y"], markers_b[1]["y"], markers_b[2]["y"] };
	double h1 = _cama_h, h2 = _camb_h;

	math::matrix<double, 9, 9> koeffs = {
		a1[0] * dcm_a["x"] + b1[0] * dcm_a["xy"] + h1 * dcm_a["xz"], 0.0, 0.0,
			-(a2[0] * dcm_b["x"] + b2[0] * dcm_b["xy"] + h1 * dcm_b["xz"]), 0.0, 0.0, -1.0, 0.0, 0.0,
		a1[0] * dcm_a["yx"] + b1[0] * dcm_a["y"] + h1 * dcm_a["yz"], 0.0, 0.0,
			-(a2[0] * dcm_b["yx"] + b2[0] * dcm_b["y"] + h1 * dcm_b["yz"]), 0.0, 0.0, 0.0, -1.0, 0.0,
		a1[0] * dcm_a["zx"] + b1[0] * dcm_a["zy"] + h1 * dcm_a["z"], 0.0, 0.0,
			-(a2[0] * dcm_b["zx"] + b2[0] * dcm_b["zy"] + h1 * dcm_b["z"]), 0.0, 0.0, 0.0, 0.0, -1.0,

		0.0, a1[1] * dcm_a["x"] + b1[1] * dcm_a["xy"] + h1 * dcm_a["xz"], 0.0,
			0.0, -(a2[1] * dcm_b["x"] + b2[1] * dcm_b["xy"] + h1 * dcm_b["xz"]), 0.0, -1.0, 0.0, 0.0,
		0.0, a1[1] * dcm_a["yx"] + b1[1] * dcm_a["y"] + h1 * dcm_a["yz"], 0.0,
			0.0, -(a2[1] * dcm_b["yx"] + b2[1] * dcm_b["y"] + h1 * dcm_b["yz"]), 0.0, 0.0, -1.0, 0.0,
		0.0, a1[1] * dcm_a["zx"] + b1[1] * dcm_a["zy"] + h1 * dcm_a["z"], 0.0,
			0.0, -(a2[1] * dcm_b["zx"] + b2[1] * dcm_b["zy"] + h1 * dcm_b["z"]), 0.0, 0.0, 0.0, -1.0,

		0.0, 0.0, a1[1] * dcm_a["x"] + b1[1] * dcm_a["xy"] + h1 * dcm_a["xz"],
			0.0, 0.0, -(a2[1] * dcm_b["x"] + b2[1] * dcm_b["xy"] + h1 * dcm_b["xz"]), -1.0, 0.0, 0.0,
		0.0, 0.0, a1[1] * dcm_a["yx"] + b1[1] * dcm_a["y"] + h1 * dcm_a["yz"],
			0.0, 0.0, -(a2[1] * dcm_b["yx"] + b2[1] * dcm_b["y"] + h1 * dcm_b["yz"]), 0.0, -1.0, 0.0,
		0.0, 0.0, a1[1] * dcm_a["zx"] + b1[1] * dcm_a["zy"] + h1 * dcm_a["z"],
			0.0, 0.0, -(a2[1] * dcm_b["zx"] + b2[1] * dcm_b["zy"] + h1 * dcm_b["z"]), 0.0, 0.0, -1.0
	};
	math::vector<double, 9> free_terms = { 0.0, 0.0, -(h1 + h2), 0.0, 0.0, -(h1 + h2), 0.0, 0.0, -(h1 + h2) };
	math::vector<double, 9> solution = math::get_slau_solution(koeffs, free_terms);
	cama_position = { 0.0, 0.0, 0.0 };
	camb_position = { solution[6], solution[7], solution[8] };
}

O3DS::proc_3d::processor_t::processor_t(const cameras_proc_data_t& cameras, const proc_2d::line_detector_t& detector,
	double dist_threshold) noexcept : _cameras(cameras), _detector(detector), _dist_threshold(dist_threshold) {}

O3DS::proc_3d::cameras_proc_data_t O3DS::proc_3d::processor_t::get_cameras_proc_data() const noexcept
{
	return _cameras;
}

void O3DS::proc_3d::processor_t::set_cameras_proc_data(const cameras_proc_data_t& cameras) noexcept
{
	_cameras = cameras;
}

O3DS::proc_2d::line_detector_t O3DS::proc_3d::processor_t::get_line_detector() const noexcept
{
	return _detector;
}

void O3DS::proc_3d::processor_t::set_line_detector(const proc_2d::line_detector_t& detector) noexcept
{
	_detector = detector;
}

void O3DS::proc_3d::processor_t::expand_point_cloud(const img_t& a, const img_t& b, std::vector<math::vec3> point_cloud) const noexcept
{
	std::list<math::vec2> points_a, points_b;
	_detector.get_line(a, points_a);
	math::mat3 dcm_a = math::dcm_12g(_cameras.cama_euler);
	double h1 = _cameras.cama_h;
	math::vec3 pos_a = _cameras.cama_position;
	_detector.get_line(b, points_b);
	math::mat3 dcm_b = math::dcm_12g(_cameras.camb_euler);
	double h2 = _cameras.camb_h;
	math::vec3 pos_b = _cameras.camb_position;
	int points_count = points_a.size();
	auto it = points_a.begin();

#pragma omp parallel for
	for (int i = 0; i < points_count; ++i)
	{
		const auto& a = *it++;
		auto l1 = [&pos_a, &h1, &dcm_a, &a](double t)
		{
			math::vec3 p = { a["x"] * t, a["y"] * t, h1 * t };
			math::vec3 o = pos_a + math::vec3(0.0, 0.0, -h1);
			return o + dcm_a * p;
		};
		for (const auto& b : points_b)
		{
			auto l2 = [&pos_b, &h2, &dcm_b, &b](double t)
			{
				math::vec3 p = { b["x"] * t, b["y"] * t, h2 * t };
				math::vec3 o = pos_b + math::vec3(0.0, 0.0, -h2);
				return o + dcm_b * p;
			};
			auto dist = [l1, l2](double t1, double t2) { return math::length(l1(t1) - l2(t2)); };
			double t1 = 1.0, t2 = 1.0;
			math::vec2 grad = math::gradient(dist, t1, t2);
			int counter = 0;
			while (math::length(grad) > _dist_threshold)
			{
				t1 -= grad["x"];
				t2 -= grad["y"];
				grad = math::gradient(dist, t1, t2);
				if (++counter >= 50) break;
			}
			if (counter < 50)
			{
#pragma omp critical
				{
					math::vec3 target_a = l1(t1);
					math::vec3 target_b = l2(t2);
					point_cloud.push_back((target_a + target_b) * 0.5);
				}
			}
		}
	}
}
