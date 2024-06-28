#include "streams.hpp"
#include "numerics.hpp"
#include "image_formats.hpp"
#include "image.hpp"

std::unordered_set<std::wstring> O3DS::img_t::_extensions = { L".bmp", L".dib" };

O3DS::img_t::img_t(const std::wstring& path) noexcept
	: _path(path)
{
	if (_extensions.find(streams::get_extension(path)) != _extensions.end())
	{
		_img = new img_formats::bmp_t;
		streams::std_binaryifstream_t stream(path);
		_load_status = _img->load_format(stream);
		stream.close();
	}
}

O3DS::img_t::img_t(size_t width, size_t height, img_formats::RGB_t fill_color) noexcept
{
	img_formats::bmp_24bit_t* img_p = new img_formats::bmp_24bit_t;
	img_p->info = new img_formats::bmp_24bit_t::common_infoblock_t(40);
	_img = img_p;
	_img->set_bitmap_size(width, height, sizeof(img_formats::RGB_t));
	_img->bitmap = new uint8_t[sizeof(img_formats::RGB_t) * height * width];
	img_formats::RGB_t* bitmap = reinterpret_cast<img_formats::RGB_t*>(_img->bitmap);
	for (int i = 0; i < width * height; ++i)
		bitmap[i] = fill_color;
}

O3DS::img_t::img_t(img_formats::format_interface_t* img) noexcept
	: _img(img) {}

O3DS::img_t::img_t(const img_t& other) noexcept
	: _load_status(other._load_status), _save_status(other._save_status), 
	_change_after_load_save(other._change_after_load_save), _path(other._path)
{
	other._img->copy(_img);
}

O3DS::img_t& O3DS::img_t::operator=(const img_t& other) noexcept
{
	if (_img) _img->clear_all_data();
	other._img->copy(_img);
	_load_status = other._load_status;
	_save_status = other._save_status;
	_change_after_load_save = other._change_after_load_save;
	_path = other._path;
	return *this;
}

O3DS::img_t::img_t(img_t&& other) noexcept
	: _img(other._img), _load_status(other._load_status), _save_status(other._save_status),
	_change_after_load_save(other._change_after_load_save), _path(other._path)
{
	other._img = nullptr;
}

O3DS::img_t& O3DS::img_t::operator=(img_t&& other) noexcept
{
	_img = other._img;
	_load_status = other._load_status;
	_save_status = other._save_status;
	_change_after_load_save = other._change_after_load_save;
	_path = other._path;
	other._img = nullptr;
	return *this;
}

O3DS::img_t::img_t() noexcept = default;

O3DS::img_t::~img_t() noexcept
{
	if (_img)
	{
		_img->clear_all_data();
		delete _img;
		_img = nullptr;
	}
}

bool O3DS::img_t::load() noexcept
{
	if (_path == L"")
		_load_status = false;
	else if (_img)
	{
		_img->clear_all_data();
		streams::std_binaryifstream_t stream(_path);
		_load_status = _img->load_format(stream);
		stream.close();
	}
	else if (!_path.empty() && _extensions.find(streams::get_extension(_path)) != _extensions.end())
	{
		_img = new img_formats::bmp_t;
		streams::std_binaryifstream_t stream(_path);
		_load_status = _img->load_format(stream);
		stream.close();
	}
	else
		_load_status = false;
	if (_load_status) _change_after_load_save = false;
	return _load_status;
}

bool O3DS::img_t::save() const noexcept
{
	if (!_img || _path == L"")
		_save_status = false;
	streams::std_binaryofstream_t stream(_path);
	if (_img->save_format(stream))
	{
		_save_status = true;
		stream.close();
	}
	else
	{
		_save_status = false;
		stream.close();
	}
	if (_save_status) _change_after_load_save = false;
	return _save_status;
}

bool O3DS::img_t::load(streams::istream_interface_t& stream) noexcept
{
	if (_img)
	{
		_img->clear_all_data();
		_load_status = _img->load_format(stream);
	}
	else
		_load_status = false;
	if (_load_status) _change_after_load_save = false;
	return _load_status;
}

bool O3DS::img_t::save(streams::ostream_interface_t& stream) const noexcept
{
	if (!_img)
		_save_status = false;
	else
		_save_status = _img->save_format(stream);
	if (_save_status) _change_after_load_save = false;
	return _save_status;
}

bool O3DS::img_t::load(const std::wstring& path) noexcept
{
	if (_img)
	{
		_img->clear_all_data();
		streams::std_binaryifstream_t stream(path);
		_load_status = _img->load_format(stream);
		stream.close();
		_path = _load_status ? path : _path;
	}
	else if (!path.empty() && _extensions.find(streams::get_extension(path)) != _extensions.end())
	{
		_img = new img_formats::bmp_t;
		streams::std_binaryifstream_t stream(path);
		_load_status = _img->load_format(stream);
		_path = _load_status ? path : _path;
		stream.close();
	}
	else
		_load_status = false;
	if (_load_status) _change_after_load_save = false;
	return _load_status;
}

bool O3DS::img_t::save(const std::wstring& path) const noexcept
{
	if (!_img)
		_save_status = false;
	streams::std_binaryofstream_t stream(path);
	if (_img->save_format(stream))
	{
		_path = path;
		_save_status = true;
		stream.close();
	}
	else
	{
		_save_status = false;
		stream.close();
	}
	if (_save_status) _change_after_load_save = false;
	return _save_status;
}

O3DS::img_formats::format_interface_t* O3DS::img_t::reset() noexcept
{
	img_formats::format_interface_t* temp = _img;
	_img = nullptr;
	_path = L"";
	_load_status = false;
	_save_status = false;
	_change_after_load_save = false;
	return temp;
}

void O3DS::img_t::set(img_formats::format_interface_t* img) noexcept
{
	if (_img)
	{
		_img->clear_all_data();
		delete _img;
		_img = img;
	}
	else
		_img = img;
	_path = L"";
	_load_status = false;
	_save_status = false;
	_change_after_load_save = true;
}

void O3DS::img_t::to_grayscale(img_formats::RGBa_t mask) noexcept
{
	if (!_img || !_img->bitmap)
		return;
	_img->to_grayscale(_img->bitmap, mask);
	_change_after_load_save = true;
}

O3DS::img_t O3DS::img_t::get_grayscale(img_formats::RGBa_t mask) const noexcept
{
	if (!_img || !_img->bitmap)
		return *this;
	O3DS::img_formats::format_interface_t* gray_img;
	_img->clone(gray_img);
	gray_img->to_grayscale(_img->bitmap, mask);
	return gray_img;
}

void O3DS::img_t::resize(size_t width, size_t height) noexcept
{
	if (!_img || !_img->bitmap)
		return;
	math::mat2 mat = {
		static_cast<double>(width) / static_cast<double>(_img->w), 0.0,
		0.0, static_cast<double>(height) / static_cast<double>(_img->h)
	};
	mat = math::inverse(mat);
	uint8_t* bitmap = new uint8_t[width * height * _img->c];
	int w = width;
	int h = height;
	int c = _img->c;
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			math::vec2 v = { j, i };
			v = mat * v;
			int x = math::round(v["x"]);
			int y = math::round(v["y"]);
			x = x < 0 ? 0 : x;
			y = y < 0 ? 0 : y;
			for (int k = 0; k < c; ++k)
				bitmap[k + j * c + i * w * c] =
					_img->bitmap[k + x * _img->c + y * _img->w * _img->c];
		}
	}
	_img->set_bitmap_size(w, h, c);
	delete[] _img->bitmap;
	_img->bitmap = bitmap;
	_change_after_load_save = true;
}

void O3DS::img_t::rotate(int deg_angle) noexcept
{
	if (!_img || !_img->bitmap)
		return;
	int mod_angle = deg_angle % 90;
	int angle = deg_angle - mod_angle;
	angle %= 360;
	int w = _img->w;
	int h = _img->h;
	uint8_t* bitmap = new uint8_t[_img->h * _img->w * _img->c];
	switch (angle)
	{
	case 90:
	{
		w = _img->h;
		h = _img->w;
#pragma omp parallel for
		for (int i = 0; i < h; ++i)
		{
			for (int j = 0; j < w; ++j)
			{
				for (int k = 0; k < _img->c; ++k)
					bitmap[k + j * _img->c + i * w * _img->c] = _img->bitmap[k + i * _img->c + j * h * _img->c];
			}
		}
		break;
	}
	case 180:
	{
#pragma omp parallel for
		for (int i = 0; i < h; ++i)
		{
			for (int j = 0; j < w; ++j)
			{
				for (int k = 0; k < _img->c; ++k)
					bitmap[k + j * _img->c + i * w * _img->c] = _img->bitmap[k + j * _img->c + (h - i - 1) * w * _img->c];
			}
		}
		break;
	}
	case 270:
	{
		w = _img->h;
		h = _img->w;
#pragma omp parallel for
		for (int i = 0; i < h; ++i)
		{
			for (int j = 0; j < w; ++j)
			{
				for (int k = 0; k < _img->c; ++k)
					bitmap[k + j * _img->c + i * w * _img->c] = _img->bitmap[k + i * _img->c + (w - j - 1) * h * _img->c];
			}
		}
		break;
	}
	}
	_img->set_bitmap_size(w, h, _img->c);
	delete[] _img->bitmap;
	_img->bitmap = bitmap;
}

O3DS::img_t O3DS::img_t::crop(math::size2 left_top, math::size2 size) const noexcept
{
	if (!_img || !_img->bitmap) return *this;
	int wt = _img->w, ht = _img->h, ct = _img->c;
	uint8_t* bitmap = _img->bitmap;
	O3DS::img_formats::format_interface_t* crop_img;
	_img->clone(crop_img);
	int w = size["x"];
	int h = size["y"];
	int c = _img->c;
	crop_img->set_bitmap_size(w, h, c);
	if (left_top["x"] + w >= _img->w || left_top["y"] + h >= _img->h)
		return crop_img;
	crop_img->bitmap = new uint8_t[w * h * c];
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < c; ++k)
				crop_img->bitmap[k + j * c + i * w * c] = _img->bitmap[k + (left_top["x"] + j) * c + (left_top["y"] + i) * _img->w * c];
		}
	}
	return crop_img;
}

void O3DS::img_t::set_path(const std::wstring& path) noexcept
{
	_path = path;
	_load_status = false;
	_save_status = false;
	_change_after_load_save = true;
}

O3DS::img_t::operator const uint8_t*() const noexcept
{
	return _img->bitmap;
}

std::wstring O3DS::img_t::get_path() const noexcept
{
	return _path;
}

uint8_t* O3DS::img_t::operator[](int index)
{
	_change_after_load_save = true;
	return _img->bitmap + index * _img->w * _img->c; 
}

const uint8_t* O3DS::img_t::operator[](int index) const
{
	return _img->bitmap + index * _img->w * _img->c;
}

bool O3DS::img_t::equal(const img_t& other) const noexcept
{
	if ((_load_status || _save_status) && (other._load_status || other._save_status)
		&& _path == other._path && !_change_after_load_save && !other._change_after_load_save) return true;
	else if (!_img || !other._img || !_img->bitmap || !other._img->bitmap || _img->w != other._img->w 
		|| _img->h != other._img->h || _img->c != other._img->c) return false;
	else
	{
		for (int i = 0; i < _img->w * _img->h * _img->c; ++i)
		{
			if (_img->bitmap[i] != other._img->bitmap[i])
				return false;
		}
		return true;
	}
}

uint8_t& O3DS::img_t::operator()(int x, int y, int z) noexcept
{
	_change_after_load_save = true;
	return _img->bitmap[z + x * _img->c + y * _img->w * _img->c];
}

uint8_t O3DS::img_t::operator()(int x, int y, int z) const noexcept
{
	return _img->bitmap[z + x * _img->c + y * _img->w * _img->c];
}

O3DS::img_formats::RGB_t& O3DS::img_t::at(int x, int y) noexcept
{
	_change_after_load_save = true;
	return *reinterpret_cast<img_formats::RGB_t*>(_img->bitmap[x + y * _img->w * _img->c]);
}

O3DS::img_formats::RGB_t O3DS::img_t::at(int x, int y) const noexcept
{
	return *reinterpret_cast<img_formats::RGB_t*>(_img->bitmap[x + y * _img->w * _img->c]);
}

size_t O3DS::img_t::width() const noexcept
{
	return _img ? _img->w : 0;
}

size_t O3DS::img_t::height() const noexcept
{
	return _img ? _img->h : 0;
}

size_t O3DS::img_t::channels() const noexcept
{
	return _img ? _img->c : 0;
}

bool O3DS::img_t::get_load_status() const noexcept
{
	return _load_status;
}

bool O3DS::img_t::get_save_status() const noexcept
{
	return _save_status;
}

const O3DS::img_formats::format_interface_t* O3DS::img_t::get_format() const noexcept
{
	return _img;
}

O3DS::img_t& O3DS::img_t::operator+=(const img_t& other) noexcept
{
	_change_after_load_save = true;
	if (!_img || !_img->bitmap || !other._img || !other._img->bitmap ||
		_img->w != other._img->w || _img->h != other._img->h || _img->c != other._img->c)
		return *this;
	int h = _img->h;
	int w = _img->w;
	int c = _img->c;
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < c; ++k)
				_img->bitmap[k + j * c + i * w * c] = std::min((int)_img->bitmap[k + j * c + i * w * c]
					+ (int)other._img->bitmap[k + j * c + i * w * c], 255);
		}
	}
	return *this;
}

O3DS::img_t& O3DS::img_t::operator-=(const img_t& other) noexcept
{
	_change_after_load_save = true;
	if (!_img || !_img->bitmap || !other._img || !other._img->bitmap ||
		_img->w != other._img->w || _img->h != other._img->h || _img->c != other._img->c)
		return *this;
	int h = _img->h;
	int w = _img->w;
	int c = _img->c;
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < c; ++k)
				_img->bitmap[k + j * c + i * w * c] = std::max((int)_img->bitmap[k + j * c + i * w * c]
					- (int)other._img->bitmap[k + j * c + i * w * c], 0);
		}
	}
	return *this;
}

O3DS::img_t O3DS::operator+(const img_t& img_a, const img_t& img_b) noexcept
{
	img_formats::format_interface_t* result;
	img_a.get_format()->clone(result);
	if (!img_a.get_format() || !img_b.get_format() ||
		img_a.width() != img_b.width() || img_a.height() != img_b.height() || img_a.channels() != img_b.channels())
		return result;
	int h = img_a.height();
	int w = img_a.width();
	int c = img_a.channels();
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < c; ++k)
				result->bitmap[k + j * c + i * w * c] = std::min((int)img_a.get_format()->bitmap[k + j * c + i * w * c]
					+ (int)img_b.get_format()->bitmap[k + j * c + i * w * c], 255);
		}
	}
	return result;
}

O3DS::img_t O3DS::operator-(const img_t& img_a, const img_t& img_b) noexcept
{
	img_formats::format_interface_t* result;
	img_a.get_format()->clone(result);
	if (!img_a.get_format() || !img_b.get_format() ||
		img_a.width() != img_b.width() || img_a.height() != img_b.height() || img_a.channels() != img_b.channels())
		return result;
	int h = img_a.height();
	int w = img_a.width();
	int c = img_a.channels();
#pragma omp parallel for
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < c; ++k)
				result->bitmap[k + j * c + i * w * c] = std::max((int)img_a.get_format()->bitmap[k + j * c + i * w * c]
					- (int)img_b.get_format()->bitmap[k + j * c + i * w * c], 0);
		}
	}
	return result;
}

bool O3DS::operator==(const img_t& img_a, const img_t& img_b) noexcept
{
	return img_a.equal(img_b);
}

bool O3DS::operator!=(const img_t& img_a, const img_t& img_b) noexcept
{
	return !(img_a == img_b);
}
