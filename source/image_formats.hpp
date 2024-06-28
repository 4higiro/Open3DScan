#ifndef _OPEN_3D_SCAN_IMG_FORMATS
#define _OPEN_3D_SCAN_IMG_FORMATS

#define DECLARE_LOAD_FUNCTION(NAME) virtual bool load_##NAME##(O3DS::streams::istream_interface_t& stream) noexcept
#define IMPLEMENT_LOAD_FUNCTION(THISTYPE, NAME) bool THISTYPE::load_##NAME##(O3DS::streams::istream_interface_t& stream) noexcept
#define OVERRIDE_LOAD_FUNCTION(NAME) bool load_##NAME##(O3DS::streams::istream_interface_t& stream) noexcept override
#define LOAD_FUNCTION_CALL(NAME) load_##NAME##(stream)
#define GET_IMG_BYTES(VAR) stream.read(reinterpret_cast<char*>(&VAR), sizeof(VAR))
#define GET_IMG_ARR(ARR, BYTES) stream.read(reinterpret_cast<char*>(ARR), BYTES)
#define DECLARE_SAVE_FUNCTION(NAME) virtual bool save_##NAME##(O3DS::streams::ostream_interface_t& stream) const noexcept
#define IMPLEMENT_SAVE_FUNCTION(THISTYPE, NAME) bool THISTYPE::save_##NAME##(O3DS::streams::ostream_interface_t& stream) const noexcept
#define OVERRIDE_SAVE_FUNCTION(NAME) bool save_##NAME##(O3DS::streams::ostream_interface_t& stream) const noexcept override
#define SAVE_FUNCTION_CALL(NAME) save_##NAME##(stream)
#define SET_IMG_BYTES(VAR) stream.write(reinterpret_cast<char*>(&VAR), sizeof(VAR));
#define SET_IMG_ARR(ARR, BYTES) stream.write(reinterpret_cast<char*>(ARR), BYTES)

using uint8_t = unsigned char;

namespace O3DS
{
	namespace streams { struct istream_interface_t; struct ostream_interface_t; }

	namespace img_formats
	{
		struct RGB_t final { uint8_t blue, green, red; RGB_t();
			RGB_t(uint8_t r, uint8_t g, uint8_t b) noexcept; };
		struct RGBa_t final { uint8_t blue, green, red, alpha; RGBa_t();
			RGBa_t(uint8_t r, uint8_t g, uint8_t b, uint8_t a) noexcept; };

		struct format_interface_t 
		{
			int w, h, c; // width, height, channel count
			uint8_t* bitmap = nullptr;

			virtual void clear_all_data() noexcept = 0;
			virtual void set_bitmap_size(int w, int h, int c) noexcept = 0;
			virtual void copy(format_interface_t*& other) const noexcept = 0;
			virtual void clone(format_interface_t*& other) const noexcept = 0;
			virtual void to_grayscale(uint8_t* from_bitmap, RGBa_t mask) noexcept = 0;

			DECLARE_LOAD_FUNCTION(format) = 0;
			DECLARE_SAVE_FUNCTION(format) = 0;
		};

		struct bmp_t : format_interface_t
		{
			struct {
				unsigned size; unsigned short reserved_no1,
					reserved_no2; unsigned offset_bits;
			} header;

			struct common_infoblock_t final
			{
				struct abstract_infoblock_t { unsigned size; };

				struct infoblock_core_t : abstract_infoblock_t
				{
					unsigned short width;
					unsigned short height;
					unsigned short plane_count;
					unsigned short bit_count;
				};

				struct infoblock_v3_t : abstract_infoblock_t
				{
					int width;
					int height;
					unsigned short plane_count;
					unsigned short bit_count;
					unsigned compression;
					unsigned size_image;
					int ppm_x, ppm_y;
					unsigned color_used;
					unsigned color_important;

					enum compression_mode_t
					{
						RGB = 0,
						RLE8 = 1,
						RLE4 = 2,
						BITFIELDS = 3,
						JPEG = 4,
						PNG = 5,
						ALPHABITFIELDS = 6
					};
				};

				struct infoblock_v4_t : infoblock_v3_t
				{
					unsigned red_mask;
					unsigned green_mask;
					unsigned blue_mask;
					unsigned alpha_mask;
					unsigned color_space_type;
					uint8_t endpoints[36];
					unsigned red_gamma;
					unsigned green_gamma;
					unsigned blue_gamma;
				};

				struct infoblock_v5_t : infoblock_v4_t
				{
					unsigned intent;
					unsigned profile_data;
					unsigned profile_size;
					unsigned reserved;
				};

				struct {
					uint8_t rno1 = 0; unsigned short rno2 = 0;
					unsigned rno3 = 0; int rno4 = 0;
				} reservs;

				abstract_infoblock_t* infoblock;

				unsigned& size;
				unsigned short& width_core;
				unsigned short& height_core;
				int& width;
				int& height;
				unsigned short& plane_count;
				unsigned short& bit_count;
				unsigned& compression;
				unsigned& size_image;
				int& ppm_x;
				int& ppm_y;
				unsigned& color_used;
				unsigned& color_important;
				unsigned& red_mask;
				unsigned& green_mask;
				unsigned& blue_mask;
				unsigned& alpha_mask;
				unsigned& color_space_type;
				uint8_t* endpoints;
				unsigned& red_gamma;
				unsigned& green_gamma;
				unsigned& blue_gamma;
				unsigned& intent;
				unsigned& profile_data;
				unsigned& profile_size;
				unsigned& reserved;

				int img_width, img_height, img_channels;

				static abstract_infoblock_t* init(unsigned size) noexcept;

				explicit common_infoblock_t(unsigned size) noexcept;
			} *info = nullptr;

			void get_most_important_color(void* palette, size_t count) const noexcept;
			void clear_all_data() noexcept override;
			void set_bitmap_size(int w, int h, int c) noexcept override;
			void copy(format_interface_t*& other) const noexcept override;
			void clone(format_interface_t*& other) const noexcept override;
			void to_grayscale(uint8_t* from_bitmap, RGBa_t mask) noexcept override;
			void from_grayscale(int channels) noexcept;

			template <typename other_bmp_format_t>
			bool save_as(streams::ostream_interface_t& stream) const noexcept
			{
				other_bmp_format_t other;
				other.header = header;
				other.info = info;
				other.bitmap = bitmap;
				return other.save_format(stream);
			}

			static constexpr int save_img_channels = 3;

			DECLARE_LOAD_FUNCTION(header);
			DECLARE_SAVE_FUNCTION(header);

			DECLARE_LOAD_FUNCTION(info);
			DECLARE_SAVE_FUNCTION(info);

			DECLARE_LOAD_FUNCTION(bitmap);
			DECLARE_SAVE_FUNCTION(bitmap);

			OVERRIDE_LOAD_FUNCTION(format);
			OVERRIDE_SAVE_FUNCTION(format);
		};

		struct bmp_mono_t : bmp_t
		{
			OVERRIDE_SAVE_FUNCTION(header);

			OVERRIDE_LOAD_FUNCTION(info);
			OVERRIDE_SAVE_FUNCTION(info);

			OVERRIDE_LOAD_FUNCTION(bitmap);
			OVERRIDE_SAVE_FUNCTION(bitmap);
		};

		struct bmp_2bit_t : bmp_t
		{
			OVERRIDE_SAVE_FUNCTION(header);

			OVERRIDE_LOAD_FUNCTION(info);
			OVERRIDE_SAVE_FUNCTION(info);

			OVERRIDE_LOAD_FUNCTION(bitmap);
			OVERRIDE_SAVE_FUNCTION(bitmap);
		};

		struct bmp_4bit_t : bmp_t
		{
			OVERRIDE_SAVE_FUNCTION(header);

			OVERRIDE_LOAD_FUNCTION(info);
			OVERRIDE_SAVE_FUNCTION(info);

			OVERRIDE_LOAD_FUNCTION(bitmap);
			OVERRIDE_SAVE_FUNCTION(bitmap);
		};

		struct bmp_8bit_t : bmp_t
		{
			OVERRIDE_SAVE_FUNCTION(header);

			OVERRIDE_LOAD_FUNCTION(info);
			OVERRIDE_SAVE_FUNCTION(info);

			OVERRIDE_LOAD_FUNCTION(bitmap);
			OVERRIDE_SAVE_FUNCTION(bitmap);
		};

		struct bmp_24bit_t : bmp_t
		{
			OVERRIDE_LOAD_FUNCTION(info);
			OVERRIDE_LOAD_FUNCTION(bitmap);
		};
	}
}

#endif
