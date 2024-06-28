#ifndef _OPEN_3D_SCAN_IMAGE
#define _OPEN_3D_SACN_IMAGE

#include "image_formats.hpp"
#include <string>
#include <list>
#include <unordered_set>

namespace O3DS
{
	namespace math { template<typename type_t, size_t dim> class vector; using size2 = vector<size_t, 2>; }

	class img_t final
	{
	private:
		mutable img_formats::format_interface_t* _img = nullptr;
		mutable std::wstring _path;
		mutable bool _load_status = false;
		mutable bool _save_status = false;
		mutable bool _change_after_load_save = false;

		static std::unordered_set<std::wstring> _extensions;

	public:
		explicit img_t(const std::wstring& path) noexcept; // Create .bmp
		img_t(size_t width, size_t height, img_formats::RGB_t fill_color = img_formats::RGB_t(0, 0, 0)) noexcept; // Create .bmp
		img_t(img_formats::format_interface_t* img) noexcept;
		img_t(const img_t& other) noexcept;
		img_t& operator=(const img_t& other) noexcept;
		img_t(img_t&& other) noexcept;
		img_t& operator=(img_t&& other) noexcept;
		img_t() noexcept;
		~img_t() noexcept;

		bool load() noexcept; // Load from get_path()
		bool save() const noexcept; // Save to get_path()
		bool load(streams::istream_interface_t& stream) noexcept;
		bool save(streams::ostream_interface_t& stream) const noexcept;
		bool load(const std::wstring& path) noexcept;
		bool save(const std::wstring& path) const noexcept;
		img_formats::format_interface_t* reset() noexcept;
		void set(img_formats::format_interface_t* img) noexcept;

		void to_grayscale(img_formats::RGBa_t mask = img_formats::RGBa_t(1, 1, 1, 1)) noexcept;
		img_t get_grayscale(img_formats::RGBa_t mask = img_formats::RGBa_t(1, 1, 1, 1)) const noexcept;
		void resize(size_t width, size_t height) noexcept;
		void rotate(int deg_angle) noexcept; // angle % 90 == 0
		img_t crop(math::size2 left_buttom, math::size2 size) const noexcept;
		void set_path(const std::wstring& path) noexcept;

		operator const uint8_t* () const noexcept;
		std::wstring get_path() const noexcept;
		uint8_t* operator[](int index);
		const uint8_t* operator[](int index) const;
		bool equal(const img_t& other) const noexcept;

		uint8_t& operator()(int x, int y, int z = 0) noexcept;
		uint8_t operator()(int x, int y, int z = 0) const noexcept;
		img_formats::RGB_t& at(int x, int y) noexcept;
		img_formats::RGB_t at(int x, int y) const noexcept;

		size_t width() const noexcept;
		size_t height() const noexcept;
		size_t channels() const noexcept;
		bool get_load_status() const noexcept;
		bool get_save_status() const noexcept;
		const img_formats::format_interface_t* get_format() const noexcept;

		img_t& operator+=(const img_t& other) noexcept;
		img_t& operator-=(const img_t& other) noexcept;
	};

	img_t operator+(const img_t& img_a, const img_t& img_b) noexcept;
	img_t operator-(const img_t& img_a, const img_t& img_b) noexcept;
	bool operator==(const img_t& img_a, const img_t& img_b) noexcept;
	bool operator!=(const img_t& img_a, const img_t& img_b) noexcept;
}

#endif
