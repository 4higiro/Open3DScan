#include "streams.hpp"
#include "image_formats.hpp"
#include <deque>
#include <algorithm>

O3DS::img_formats::RGB_t::RGB_t() = default;

O3DS::img_formats::RGB_t::RGB_t(uint8_t r, uint8_t g, uint8_t b) noexcept
	: red(r), green(g), blue(b) {}

O3DS::img_formats::RGBa_t::RGBa_t() = default;

O3DS::img_formats::RGBa_t::RGBa_t(uint8_t r, uint8_t g, uint8_t b, uint8_t a) noexcept
	: red(r), green(g), blue(b), alpha(a) {}

#define INIT_FIELD(FIELD, RNO) (size == sizeof(infoblock_core_t) ? reinterpret_cast<infoblock_core_t*>(infoblock)->FIELD : \
						(size == sizeof(infoblock_v3_t) ? reinterpret_cast<infoblock_v3_t*>(infoblock)->FIELD : \
						(size == sizeof(infoblock_v4_t) ? reinterpret_cast<infoblock_v4_t*>(infoblock)->FIELD : \
						(size == sizeof(infoblock_v5_t) ? reinterpret_cast<infoblock_v5_t*>(infoblock)->FIELD : reservs.RNO))))
#define INIT_UNIQ_FIELD(TYPE, FIELD, RNO) (size == sizeof(TYPE) ? reinterpret_cast<TYPE*>(infoblock)->FIELD : reservs.RNO)
#define INIT_OF2_FIELD(TYPEA, TYPEB, FIELD, RNO) (size == sizeof(TYPEA) ? reinterpret_cast<TYPEA*>(infoblock)->FIELD : \
						(size == sizeof(TYPEB) ? reinterpret_cast<TYPEB*>(infoblock)->FIELD : reservs.RNO))
#define INIT_OF3_FIELD(TYPEA, TYPEB, TYPEC, FIELD, RNO) (size == sizeof(TYPEA) ? reinterpret_cast<TYPEA*>(infoblock)->FIELD : \
						(size == sizeof(TYPEB) ? reinterpret_cast<TYPEB*>(infoblock)->FIELD : \
						(size == sizeof(TYPEC) ? reinterpret_cast<TYPEC*>(infoblock)->FIELD : reservs.RNO)))

O3DS::img_formats::bmp_t::common_infoblock_t::abstract_infoblock_t* 
	O3DS::img_formats::bmp_t::common_infoblock_t::init(unsigned size) noexcept
{
	switch (size)
		{
		case sizeof(infoblock_core_t) :
			return new infoblock_core_t;
		case sizeof(infoblock_v3_t) :
			return new infoblock_v3_t;
		case sizeof(infoblock_v4_t) :
			return new infoblock_v4_t;
		case sizeof(infoblock_v5_t) :
			return new infoblock_v5_t;
		default:
			return new abstract_infoblock_t;
	}
}

O3DS::img_formats::bmp_t::common_infoblock_t::common_infoblock_t(unsigned size) noexcept
	: infoblock(init(size)), size(infoblock->size),
	width_core(INIT_UNIQ_FIELD(infoblock_core_t, width, rno2)),
	height_core(INIT_UNIQ_FIELD(infoblock_core_t, height, rno2)),
	width(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, width, rno4)),
	height(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, height, rno4)),
	plane_count(INIT_FIELD(plane_count, rno2)), bit_count(INIT_FIELD(bit_count, rno2)),
	compression(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, compression, rno3)),
	size_image(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, size_image, rno3)),
	ppm_x(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, ppm_x, rno4)),
	ppm_y(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, ppm_y, rno4)),
	color_used(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, color_used, rno3)),
	color_important(INIT_OF3_FIELD(infoblock_v3_t, infoblock_v4_t, infoblock_v5_t, color_important, rno3)),
	red_mask(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, red_mask, rno3)),
	green_mask(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, green_mask, rno3)),
	blue_mask(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, blue_mask, rno3)),
	alpha_mask(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, alpha_mask, rno3)),
	color_space_type(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, color_space_type, rno3)),
	red_gamma(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, red_gamma, rno3)),
	green_gamma(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, green_gamma, rno3)),
	blue_gamma(INIT_OF2_FIELD(infoblock_v4_t, infoblock_v5_t, blue_gamma, rno3)),
	intent(INIT_UNIQ_FIELD(infoblock_v5_t, intent, rno3)),
	profile_data(INIT_UNIQ_FIELD(infoblock_v5_t, profile_data, rno3)),
	profile_size(INIT_UNIQ_FIELD(infoblock_v5_t, profile_size, rno3)),
	reserved(INIT_UNIQ_FIELD(infoblock_v5_t, reserved, rno3))
{
	endpoints = (size == sizeof(infoblock_v4_t) ? reinterpret_cast<infoblock_v4_t*>(infoblock)->endpoints
		: (size == sizeof(infoblock_v5_t) ? reinterpret_cast<infoblock_v5_t*>(infoblock)->endpoints : nullptr));
}

void O3DS::img_formats::bmp_t::get_most_important_color(void* palette, size_t count) const noexcept
{
	if (!palette || count == 0) return;
	const int h = info->img_height;
	const int w = info->img_width;
	const int n = info->img_channels;
	switch (n)
	{
	case 3:
	{
#define IMPORTANT_BITS(TYPE) TYPE* p = reinterpret_cast<TYPE*>(palette);									\
		TYPE* b = reinterpret_cast<TYPE*>(bitmap);															\
		size_t* c = new size_t[w * h];																		\
		for (int i = 0; i < h * w; ++i)																		\
		{																									\
			size_t count = 1;																				\
			TYPE color = b[i];																				\
			for (int j = i - 1; j >= 0; ++j)																\
			{																								\
				TYPE current = b[j];																		\
				if (current.red == color.red && current.green == color.green && current.blue == color.blue)	\
				{																							\
					count += c[j];																			\
					break;																					\
				}																							\
			}																								\
			c[i] = count;																					\
		}																									\
		std::deque<size_t> max_indices;																		\
		size_t max_index = 0;																				\
		for (int i = 0; i < h * w; ++i)																		\
		{																									\
			if (c[i] > c[max_index] && find_if(max_indices.begin(), max_indices.end(),						\
				[i, b](size_t index)																		\
				{																							\
					return b[i].red == b[index].red															\
						&& b[i].green == b[index].green														\
						&& b[i].blue == b[index].blue;														\
				}) == max_indices.end())																	\
			{																								\
				max_index = i;																				\
				max_indices.push_back(max_index);															\
				if (max_indices.size() > count)																\
					max_indices.pop_front();																\
			}																								\
		}																									\
		for (int i = 0; i < max_indices.size(); ++i)														\
			p[i] = b[max_indices[i]];																		
		IMPORTANT_BITS(RGB_t)
		for (int i = max_indices.size(); i < count; ++i)													
			p[i] = RGB_t(0, 0, 0);
	}
	case 4:
	{
		IMPORTANT_BITS(RGBa_t)
		for (int i = max_indices.size(); i < count; ++i)													
			p[i] = RGBa_t(0, 0, 0, 0);
	}
	default:
		return;
	}
}

void O3DS::img_formats::bmp_t::clear_all_data() noexcept
{
	if (bitmap) delete[] bitmap;
	header.size = 0;
	header.reserved_no1 = 0;
	header.reserved_no2 = 0;
	header.offset_bits = 0;
	if (info) delete info;
}

void O3DS::img_formats::bmp_t::set_bitmap_size(int w, int h, int c) noexcept
{
	this->w = w; this->h = h; this->c = c;
	info->img_width = w;
	info->img_height = h;
	info->img_channels = c;
}

void O3DS::img_formats::bmp_t::copy(format_interface_t*& other) const noexcept
{
	other = new img_formats::bmp_t;
	auto copy_to_ptr = dynamic_cast<bmp_t*>(other);
	if (copy_to_ptr && info && bitmap)
	{
		copy_to_ptr->header = header;
		copy_to_ptr->info = new common_infoblock_t(info->size);
		copy_to_ptr->info->width_core = info->width_core;
		copy_to_ptr->info->height_core = info->height_core;
		copy_to_ptr->info->width = info->width;
		copy_to_ptr->info->height = info->height;
		copy_to_ptr->info->plane_count = info->plane_count;
		copy_to_ptr->info->bit_count = info->bit_count;
		copy_to_ptr->info->compression = info->compression;
		copy_to_ptr->info->size_image = info->size_image;
		copy_to_ptr->info->ppm_x = info->ppm_x;
		copy_to_ptr->info->ppm_y = info->ppm_y;
		copy_to_ptr->info->color_used = info->color_used;
		copy_to_ptr->info->color_important = info->color_important;
		copy_to_ptr->info->red_mask = info->red_mask;
		copy_to_ptr->info->green_mask = info->green_mask;
		copy_to_ptr->info->blue_mask = info->blue_mask;
		copy_to_ptr->info->alpha_mask = info->alpha_mask;
		copy_to_ptr->info->color_space_type = info->color_space_type;
		copy_to_ptr->info->red_gamma = info->red_gamma;
		copy_to_ptr->info->green_gamma = info->green_gamma;
		copy_to_ptr->info->blue_gamma = info->blue_gamma;
		copy_to_ptr->info->intent = info->intent;
		copy_to_ptr->info->profile_data = info->profile_data;
		copy_to_ptr->info->profile_size = info->profile_size;
		copy_to_ptr->info->reserved = info->reserved;
		copy_to_ptr->set_bitmap_size(w, h, c);
		copy_to_ptr->bitmap = new uint8_t[w * h * c];
		int size = w * h * c;
		for (int i = 0; i < size; ++i)
			copy_to_ptr->bitmap[i] = bitmap[i];
	}
}

void O3DS::img_formats::bmp_t::clone(format_interface_t*& other) const noexcept
{
	other = new img_formats::bmp_t;
	auto clone_to_ptr = dynamic_cast<bmp_t*>(other);
	if (clone_to_ptr && info && bitmap)
	{
		clone_to_ptr->header = header;
		clone_to_ptr->info = new common_infoblock_t(info->size);
		clone_to_ptr->info->width_core = info->width_core;
		clone_to_ptr->info->height_core = info->height_core;
		clone_to_ptr->info->width = info->width;
		clone_to_ptr->info->height = info->height;
		clone_to_ptr->info->plane_count = info->plane_count;
		clone_to_ptr->info->bit_count = info->bit_count;
		clone_to_ptr->info->compression = info->compression;
		clone_to_ptr->info->size_image = info->size_image;
		clone_to_ptr->info->ppm_x = info->ppm_x;
		clone_to_ptr->info->ppm_y = info->ppm_y;
		clone_to_ptr->info->color_used = info->color_used;
		clone_to_ptr->info->color_important = info->color_important;
		clone_to_ptr->info->red_mask = info->red_mask;
		clone_to_ptr->info->green_mask = info->green_mask;
		clone_to_ptr->info->blue_mask = info->blue_mask;
		clone_to_ptr->info->alpha_mask = info->alpha_mask;
		clone_to_ptr->info->color_space_type = info->color_space_type;
		clone_to_ptr->info->red_gamma = info->red_gamma;
		clone_to_ptr->info->green_gamma = info->green_gamma;
		clone_to_ptr->info->blue_gamma = info->blue_gamma;
		clone_to_ptr->info->intent = info->intent;
		clone_to_ptr->info->profile_data = info->profile_data;
		clone_to_ptr->info->profile_size = info->profile_size;
		clone_to_ptr->info->reserved = info->reserved;
		clone_to_ptr->set_bitmap_size(w, h, c);
		clone_to_ptr->bitmap = new uint8_t[w * h * c];
	}
}

void O3DS::img_formats::bmp_t::to_grayscale(uint8_t* from_bitmap, RGBa_t mask) noexcept
{
	if (!info || !bitmap || c == 1) return;
	uint8_t* gray_bitmap = new uint8_t[w * h];
	uint8_t* mask_ptr = reinterpret_cast<uint8_t*>(&mask);
	int sum = (bool)mask.blue + (bool)mask.green + (bool)mask.red;
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			short color = 0;
			for (int k = 0; k < c; ++k)
				color += mask_ptr[k] > 0 ? from_bitmap[k + j * c + i * w * c] : 0;
			gray_bitmap[j + i * w] = color / sum;
		}
	}
	delete[] bitmap;
	bitmap = gray_bitmap;
	set_bitmap_size(w, h, 1);
}

void O3DS::img_formats::bmp_t::from_grayscale(int channels) noexcept
{
	if (!info || !bitmap) return;
	uint8_t* color_bitmap = new uint8_t[w * h * channels];
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			for (int k = 0; k < channels; ++k)
			{ 
				if (k == 4) bitmap[j + i * w] = 0;
				color_bitmap[k + j * channels + i * w * channels] = bitmap[j + i * w];
			}
		}
	}
	delete bitmap;
	bitmap = color_bitmap;
	set_bitmap_size(w, h, channels);
#define FROM_GRAYSCALE(C, TYPE) if (c == 1) const_cast<TYPE*>(this)->from_grayscale(C);
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_t, header)
{
	unsigned short img_type;
	GET_IMG_BYTES(img_type);
	if (img_type != 19778)
		return false;
	GET_IMG_BYTES(header);
	if (header.reserved_no1 != 0 || header.reserved_no2 != 0)
		return false;
	return true;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_t, header)
{
#define SAVE_HEADER_BITS(BITS) if (!bitmap || !info || !stream.is_open())											\
		return false;																								\
	unsigned short img_type = 19778;																				\
	SET_IMG_BYTES(img_type);																						\
	int row_bytes = info->img_width * save_img_channels;															\
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;															\
	unsigned size = sizeof(unsigned short) + sizeof(header) + sizeof(common_infoblock_t::infoblock_v3_t)			\
		+ info->img_height * (row_bytes + need);																	\
	SET_IMG_BYTES(size);																							\
	unsigned short reserved_no1 = 0, reserved_no2 = 0;																\
	SET_IMG_BYTES(reserved_no1);																					\
	SET_IMG_BYTES(reserved_no2);																					\
	unsigned offset_bits = sizeof(unsigned short) + sizeof(header)													\
		+ sizeof(common_infoblock_t::infoblock_v3_t) + BITS * sizeof(RGBa_t);										\
	SET_IMG_BYTES(offset_bits);
	SAVE_HEADER_BITS(0)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_t, info)
{
#define LOAD_INFO_BITS if (header.size == 0)									\
		return false;															\
	unsigned size;																\
	GET_IMG_BYTES(size);														\
	info = new common_infoblock_t(size);										\
	stream.seekg(sizeof(unsigned short) + sizeof(header), std::ios_base::beg);	\
	GET_IMG_ARR(info->infoblock, size);											\
	if (size == sizeof(common_infoblock_t::infoblock_core_t))					\
		info->img_width = info->width_core,										\
		info->img_height = info->height_core;									\
	else																		\
		info->img_width = abs(info->width),										\
		info->img_height = abs(info->height);									\
	w = info->img_width;														\
	h = info->img_height;														\
	c = info->img_channels;
	LOAD_INFO_BITS
	return info->bit_count > 0 && info->compression == 0;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_t, info)
{
	FROM_GRAYSCALE(3, bmp_t)
	int row_bytes = info->img_width * save_img_channels;
#define SAVE_INFO_BITS(BIT_COUNT) int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;		\
	if (!bitmap || !info || !stream.is_open())													\
		return false;																			\
	common_infoblock_t::infoblock_v3_t save_infoblock;											\
	save_infoblock.size = sizeof(common_infoblock_t::infoblock_v3_t);							\
	save_infoblock.width = info->img_width;														\
	save_infoblock.height = info->img_height;													\
	save_infoblock.plane_count = 1;																\
	save_infoblock.bit_count = BIT_COUNT;														\
	save_infoblock.compression = common_infoblock_t::infoblock_v3_t::compression_mode_t::RGB;	\
	save_infoblock.size_image = (row_bytes + need) * info->img_height;							\
	save_infoblock.ppm_x = save_infoblock.ppm_y = 0;											\
	save_infoblock.color_used = 0;																\
	save_infoblock.color_important = 0;															\
	SET_IMG_BYTES(save_infoblock);																
	SAVE_INFO_BITS(24)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_t, bitmap)
{
	if (!info) return false;
	unsigned red_mask, green_mask, blue_mask, alpha_mask = 0;
	bool mask_used;
	if (mask_used = info->size == sizeof(common_infoblock_t::infoblock_v3_t) ? (info->bit_count % 16 == 0
		&& info->compression % 3 == 0 && info->bit_count != 64) : false)
	{
		GET_IMG_BYTES(red_mask);
		GET_IMG_BYTES(green_mask);
		GET_IMG_BYTES(blue_mask);
		if (info->compression == 6) GET_IMG_BYTES(alpha_mask);
	}
	void* palette = nullptr; int palette_size = 0; unsigned rgb_size = 0;
	if (info->size == sizeof(common_infoblock_t::infoblock_core_t) && info->bit_count <= 8)
	{
		palette_size = std::pow(2.0, info->bit_count);
		palette = new uint8_t[palette_size * sizeof(RGB_t)];
		GET_IMG_ARR(palette, palette_size * sizeof(RGB_t));
		rgb_size = sizeof(RGB_t);
	}
	else if (info->size >= sizeof(common_infoblock_t::infoblock_v3_t) && info->bit_count <= 8)
	{
		palette_size = info->color_used ? info->color_used : std::pow(2.0, info->bit_count);
		palette = new uint8_t[palette_size * sizeof(RGBa_t)];
		GET_IMG_ARR(palette, palette_size * sizeof(RGBa_t));
		rgb_size = sizeof(RGBa_t);
	}
	else if (info->size >= sizeof(common_infoblock_t::infoblock_v3_t) && info->bit_count > 8 && info->color_used)
	{
		palette_size = info->color_used;
		palette = new uint8_t[palette_size * sizeof(RGBa_t)];
		GET_IMG_ARR(palette, palette_size * sizeof(RGBa_t));
		rgb_size = sizeof(RGBa_t);
	}
	info->img_channels = rgb_size ? rgb_size : sizeof(RGB_t);
	c = info->img_channels;
	switch (info->bit_count)
	{
	case 1: case 2: case 4:
	{
		if (palette_size > 0 && rgb_size == 4)
		{
			stream.seekg(header.offset_bits, std::ios_base::beg);
			RGBa_t* p = reinterpret_cast<RGBa_t*>(palette);
			RGBa_t* arr = reinterpret_cast<RGBa_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGBa_t)]);
#define BIT_1_2_4_LOAD  const int divider = 8 / info->bit_count;													\
						const int row_bytes = info->img_width / divider + (info->img_width % divider == 0 ? 0 : 1);	\
						const int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;								\
						for (int i = 0; i < info->img_height; ++i)													\
						{																							\
							int counter = 0;																		\
							for (int j = 0; j < row_bytes; ++j)														\
							{																						\
								uint8_t byte;																		\
								GET_IMG_BYTES(byte);																\
								bool bits[8];																		\
								for (int k = 0; k < 8; ++k)															\
								{																					\
									bits[k] = byte % 2;																\
									byte >>= 1;																		\
								}																					\
								for (int k = 0; k < 8; k += info->bit_count)										\
								{																					\
									uint8_t index = bits[7 - k];													\
									for (int l = k + 1; l < k + info->bit_count; ++l)								\
									{																				\
										index <<= 1;																\
										index += bits[7 - l];														\
									}																				\
									if (counter < info->img_width)													\
									{																				\
										arr[counter + (info->img_height - 1 - i) * info->img_width] = p[index];		\
										++counter;																	\
									}																				\
									else																			\
										break;																		\
								}																					\
							}																						\
							uint8_t buff[4];																		\
							GET_IMG_ARR(buff, need);																\
						}																							\
						bitmap = reinterpret_cast<uint8_t*>(arr);													
			BIT_1_2_4_LOAD
		}
		else if (palette_size > 0 && rgb_size == 3)
		{
			stream.seekg(header.offset_bits, std::ios_base::beg);
			RGB_t* p = reinterpret_cast<RGB_t*>(palette);
			RGB_t* arr = reinterpret_cast<RGB_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGB_t)]);
			BIT_1_2_4_LOAD
		}
		else
			bitmap = new uint8_t[info->img_width * info->img_height * info->img_channels]();
		break;
	}
	case 8:
	{
		if (palette_size > 0 && rgb_size == 4)
		{
			stream.seekg(header.offset_bits, std::ios_base::beg);
			RGBa_t* p = reinterpret_cast<RGBa_t*>(palette);
			RGBa_t* arr = reinterpret_cast<RGBa_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGBa_t)]);
#define BIT_8_LOAD		const int need = info->img_width % 4 == 0 ? 0 : 4 - info->img_width % 4;	\
						for (int i = 0; i < info->img_height; ++i)									\
						{																			\
							for (int j = 0; j < info->img_width; ++j)								\
							{																		\
								uint8_t index;														\
								GET_IMG_BYTES(index);												\
								arr[j + (info->img_height - 1 - i) * info->img_width] = p[index];	\
							}																		\
							uint8_t buff[4];														\
							GET_IMG_ARR(buff, need);												\
						}																			\
						bitmap = reinterpret_cast<uint8_t*>(arr);	
			BIT_8_LOAD
		}
		else if (palette_size > 0 && rgb_size == 3)
		{
			stream.seekg(header.offset_bits, std::ios_base::beg);
			RGB_t* p = reinterpret_cast<RGB_t*>(palette);
			RGB_t* arr = reinterpret_cast<RGB_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGB_t)]);
			BIT_8_LOAD
		}
		else
			bitmap = new uint8_t[info->img_width * info->img_height * info->img_channels]();
		break;
	}
	case 24:
	{
		RGB_t* arr = new RGB_t[info->img_width * info->img_height];
		const int need = info->img_width * sizeof(RGB_t) % 4 == 0 ? 0 : 4 - (info->img_width * sizeof(RGB_t) % 4);
		const int bias = info->img_width * sizeof(RGB_t) + need;
		for (int i = 0; i < info->img_height; ++i)
		{
			uint8_t buff[4];
			stream.seekg(header.offset_bits + i * bias, std::ios_base::beg);
			GET_IMG_ARR(arr + (info->img_height - 1 - i) * info->img_width, info->img_width * sizeof(RGB_t));
			GET_IMG_ARR(buff, need);
		}
		bitmap = reinterpret_cast<uint8_t*>(arr);
		break;
	}
	}
	return bitmap;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_t, bitmap)
{
	FROM_GRAYSCALE(3, bmp_t)
	int row_bytes = info->img_width * save_img_channels;
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;
	static uint8_t buff[4] = { 0, 0, 0, 0 };
	if (!bitmap || !info || !stream.is_open())
		return false;
	else if (info->img_channels == save_img_channels)
	{
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * row_bytes;
		for (int i = 0; i < info->img_height; ++i)
		{
			SET_IMG_ARR(save_bmp, row_bytes);
			SET_IMG_ARR(buff, need);
			save_bmp -= row_bytes;
		}
		return true;
	}
	else if (info->img_channels > save_img_channels)
	{
		int real_row_bytes = info->img_width * info->img_channels;
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * real_row_bytes;
		for (int i = 0; i < info->img_height; ++i)
		{
			uint8_t* save_row = save_bmp;
			for (int j = 0; j < info->img_width; ++j)
			{
				SET_IMG_ARR(save_row, save_img_channels);
				save_row += info->img_channels;
			}
			SET_IMG_ARR(buff, need);
			save_bmp -= real_row_bytes;
		}
		return true;
	}
	else
		return false;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_t, format)
{
	if (!LOAD_FUNCTION_CALL(header))
		return false;
	if (!LOAD_FUNCTION_CALL(info))
		return false;
	if (!LOAD_FUNCTION_CALL(bitmap))
		return false;
	return true;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_t, format)
{
	if (!SAVE_FUNCTION_CALL(header))
		return false;
	if (!SAVE_FUNCTION_CALL(info))
		return false;
	if (!SAVE_FUNCTION_CALL(bitmap))
		return false;
	return true;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_mono_t, header)
{
	SAVE_HEADER_BITS(2)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_mono_t, info)
{
	LOAD_INFO_BITS
	return info->bit_count == 1 && info->compression == 0;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_mono_t, info)
{
	FROM_GRAYSCALE(4, bmp_mono_t)
	int row_bytes = info->img_width / 8 + info->img_width % 8 == 0 ? 0 : 1;
	SAVE_INFO_BITS(1)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_mono_t, bitmap)
{
	if (!info || info->bit_count != 1)
		return false;
#define LOAD_BITMAP_BITS(BITS) void* palette = nullptr; int palette_size = 0; unsigned rgb_size = 0;				\
	if (info->size == sizeof(common_infoblock_t::infoblock_core_t))													\
	{																												\
		palette_size = std::pow(2.0, info->bit_count);																\
		palette = new uint8_t[palette_size * sizeof(RGB_t)];														\
		GET_IMG_ARR(palette, palette_size * sizeof(RGB_t));															\
		rgb_size = sizeof(RGB_t);																					\
	}																												\
	else if (info->size >= sizeof(common_infoblock_t::infoblock_v3_t))												\
	{																												\
		palette_size = info->color_used ? info->color_used : std::pow(2.0, info->bit_count);						\
		palette = new uint8_t[palette_size * sizeof(RGBa_t)];														\
		GET_IMG_ARR(palette, palette_size * sizeof(RGBa_t));														\
		rgb_size = sizeof(RGBa_t);																					\
	}																												\
	else if (info->size >= sizeof(common_infoblock_t::infoblock_v3_t) && info->color_used)							\
	{																												\
		palette_size = info->color_used;																			\
		palette = new uint8_t[palette_size * sizeof(RGBa_t)];														\
		GET_IMG_ARR(palette, palette_size * sizeof(RGBa_t));														\
		rgb_size = sizeof(RGBa_t);																					\
	}																												\
	info->img_channels = rgb_size ? rgb_size : sizeof(RGB_t);														\
	c = info->img_channels;																							\
	if (palette_size > 0 && rgb_size == 4)																			\
	{																												\
		stream.seekg(header.offset_bits, std::ios_base::beg);														\
		RGBa_t* p = reinterpret_cast<RGBa_t*>(palette);																\
		RGBa_t* arr = reinterpret_cast<RGBa_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGBa_t)]);	\
		BITS																										\
	}																												\
	else if (palette_size > 0 && rgb_size == 3)																		\
	{																												\
		stream.seekg(header.offset_bits, std::ios_base::beg);														\
		RGB_t* p = reinterpret_cast<RGB_t*>(palette);																\
		RGB_t* arr = reinterpret_cast<RGB_t*>(new uint8_t[info->img_width * info->img_height * sizeof(RGB_t)]);		\
		BITS																										\
	}																												\
	else																											\
		bitmap = new uint8_t[info->img_width * info->img_height * info->img_channels]();
	LOAD_BITMAP_BITS(BIT_1_2_4_LOAD)
	return bitmap;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_mono_t, bitmap)
{
	FROM_GRAYSCALE(4, bmp_mono_t)
	int row_bytes = info->img_width / 8 + info->img_width % 8 == 0 ? 0 : 1;
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;
	static uint8_t buff[4] = { 0, 0, 0, 0 };
	if (!bitmap || !info || !stream.is_open())
		return false;
	else if (info->img_channels == save_img_channels)
	{
#define SAVE_1_BIT(TYPE) RGBa_t palette[2] = { RGBa_t(0, 0, 0, 0), RGBa_t(255, 255, 255, 0) };				\
		SET_IMG_ARR(palette, 2 * sizeof(RGBa_t));															\
		uint8_t r1 = palette[0].red;																		\
		uint8_t g1 = palette[0].green;																		\
		uint8_t b1 = palette[0].blue;																		\
		uint8_t r2 = palette[1].red;																		\
		uint8_t g2 = palette[1].green;																		\
		uint8_t b2 = palette[1].blue;																		\
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * (info->img_width * info->img_channels);		\
		for (int i = 0; i < info->img_height; ++i)															\
		{																									\
			TYPE* p = reinterpret_cast<TYPE*>(save_bmp);													\
			for (int j = 0; j < info->img_width;)															\
			{																								\
				int k0 = -1;																				\
				uint8_t index = 0b00'000'000;																\
				for (int k = 0; k < 8; ++k)																	\
				{																							\
					if (j >= info->img_width)																\
					{																						\
						k0 = k - 1;																			\
						break;																				\
					}																						\
					uint8_t r0 = p[j].red;																	\
					uint8_t g0 = p[j].green;																\
					uint8_t b0 = p[j].blue;																	\
					double l1 = sqrt((r1 - r0) * (r1 - r0) + (g1 - g0) * (g1 - g0) + (b1 - b0) * (b1 - b0));\
					double l2 = sqrt((r2 - r0) * (r2 - r0) + (g2 - g0) * (g2 - g0) + (b2 - b0) * (b2 - b0));\
					if (l1 > l2) index += pow(2, 7 - k);													\
					++j;																					\
				}																							\
				SET_IMG_BYTES(index);																		\
			}																								\
			SET_IMG_ARR(buff, need);																		\
			save_bmp -= info->img_width * info->img_channels;												\
		}
		SAVE_1_BIT(RGB_t)
		return true;
	}
	else if (info->img_channels > save_img_channels)
	{
		SAVE_1_BIT(RGBa_t)
		return true;
	}
	else
		return false;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_2bit_t, header)
{
	SAVE_HEADER_BITS(4)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_2bit_t, info)
{
	LOAD_INFO_BITS
	return info->bit_count == 2 && info->compression == 0;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_2bit_t, info)
{
	FROM_GRAYSCALE(4, bmp_2bit_t)
	int row_bytes = info->img_width / 4 + info->img_width % 8 == 0 ? 0 : 1;
	SAVE_INFO_BITS(2)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_2bit_t, bitmap)
{
	if (!info || info->bit_count != 2) return false;
	LOAD_BITMAP_BITS(BIT_1_2_4_LOAD)
	return bitmap;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_2bit_t, bitmap)
{
	FROM_GRAYSCALE(4, bmp_2bit_t)
	int row_bytes = info->img_width / 4 + info->img_width % 8 == 0 ? 0 : 1;
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;
	static uint8_t buff[4] = { 0, 0, 0, 0 };
	if (!bitmap || !info || !stream.is_open())
		return false;
	else if (info->img_channels == save_img_channels)
	{
#define SAVE_2_BIT(TYPE) RGBa_t palette[4];																\
		get_most_important_color(palette, 4);															\
		SET_IMG_ARR(palette, 4 * sizeof(RGBa_t));														\
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * (info->img_width * info->img_channels);	\
		for (int i = 0; i < info->img_height; ++i)														\
		{																								\
			TYPE* p = reinterpret_cast<TYPE*>(save_bmp);												\
			for (int j = 0;;)																			\
			{																							\
				uint8_t index = 0, min_index = 0;														\
				double min_len = 255.0 * 255.0 * 255.0;													\
				for (int k = 0; k < 4; ++k)																\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				++j;																					\
				index = min_index;																		\
				if (j >= info->img_width)																\
				{																						\
					index <<= 6;																		\
					SET_IMG_BYTES(index);																\
					break;																				\
				}																						\
				min_index = 0; min_len = 255.0 * 255.0 * 255.0;											\
				for (int k = 0; k < 4; ++k)																\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				index <<= 2;																			\
				index += min_index;																		\
				++j;																					\
				if (j >= info->img_width)																\
				{																						\
					index <<= 4;																		\
					SET_IMG_BYTES(index);																\
					break;																				\
				}																						\
				min_index = 0; min_len = 255.0 * 255.0 * 255.0;											\
				for (int k = 0; k < 4; ++k)																\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				index <<= 2;																			\
				index += min_index;																		\
				++j;																					\
				if (j >= info->img_width)																\
				{																						\
					index <<= 2;																		\
					SET_IMG_BYTES(index);																\
					break;																				\
				}																						\
				min_index = 0; min_len = 255.0 * 255.0 * 255.0;											\
				for (int k = 0; k < 4; ++k)																\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				index <<= 2;																			\
				index += min_index;																		\
				++j;																					\
				SET_IMG_BYTES(index);																	\
				if (j >= info->img_width)																\
					break;																				\
			}																							\
			SET_IMG_ARR(buff, need);																	\
			save_bmp -= info->img_width * info->img_channels;											\
		}
		SAVE_2_BIT(RGB_t)
			return true;
	}
	else if (info->img_channels > save_img_channels)
	{
		SAVE_2_BIT(RGBa_t)
			return true;
	}
	else
		return false;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_4bit_t, header)
{
	SAVE_HEADER_BITS(16)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_4bit_t, info)
{
	LOAD_INFO_BITS
	return info->bit_count == 4 && info->compression == 0;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_4bit_t, info)
{
	FROM_GRAYSCALE(4, bmp_4bit_t)
	int row_bytes = info->img_width / 2 + info->img_width % 8 == 0 ? 0 : 1;
	SAVE_INFO_BITS(4)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_4bit_t, bitmap)
{
	if (!info || info->bit_count != 4) return false;
	LOAD_BITMAP_BITS(BIT_1_2_4_LOAD)
	return bitmap;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_4bit_t, bitmap)
{
	FROM_GRAYSCALE(4, bmp_4bit_t)
	int row_bytes = info->img_width / 2 + info->img_width % 8 == 0 ? 0 : 1;
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;
	static uint8_t buff[4] = { 0, 0, 0, 0 };
	if (!bitmap || !info || !stream.is_open())
		return false;
	else if (info->img_channels == save_img_channels)
	{
#define SAVE_4_BIT(TYPE) RGBa_t palette[16];															\
		get_most_important_color(palette, 16);															\
		SET_IMG_ARR(palette, 16 * sizeof(RGBa_t));														\
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * (info->img_width * info->img_channels);	\
		for (int i = 0; i < info->img_height; ++i)														\
		{																								\
			TYPE* p = reinterpret_cast<TYPE*>(save_bmp);												\
			for (int j = 0;;)																			\
			{																							\
				uint8_t index = 0, min_index = 0;														\
				double min_len = 255.0 * 255.0 * 255.0;													\
				for (int k = 0; k < 16; ++k)															\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				++j;																					\
				index = min_index;																		\
				if (j >= info->img_width)																\
				{																						\
					index <<= 4;																		\
					SET_IMG_BYTES(index);																\
					break;																				\
				}																						\
				min_index = 0; min_len = 255.0 * 255.0 * 255.0;											\
				for (int k = 0; k < 16; ++k)															\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				index <<= 4;																			\
				index += min_index;																		\
				++j;																					\
				SET_IMG_BYTES(index);																	\
				if (j >= info->img_width)																\
					break;																				\
			}																							\
			SET_IMG_ARR(buff, need);																	\
			save_bmp -= info->img_width * info->img_channels;											\
		}
		SAVE_4_BIT(RGB_t)
		return true;
	}
	else if (info->img_channels > save_img_channels)
	{
		SAVE_4_BIT(RGBa_t)
		return true;
	}
	else
		return false;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_8bit_t, header)
{
	SAVE_HEADER_BITS(256)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_8bit_t, info)
{
	LOAD_INFO_BITS
	return info->bit_count == 8 && info->compression == 0;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_8bit_t, info)
{
	FROM_GRAYSCALE(4, bmp_8bit_t)
	int row_bytes = info->img_width;
	SAVE_INFO_BITS(8)
	return true;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_8bit_t, bitmap)
{
	if (!info || info->bit_count != 8) return false;
	LOAD_BITMAP_BITS(BIT_8_LOAD)
	return bitmap;
}

IMPLEMENT_SAVE_FUNCTION(O3DS::img_formats::bmp_8bit_t, bitmap)
{
	FROM_GRAYSCALE(4, bmp_8bit_t)
	int row_bytes = info->img_width;
	int need = row_bytes % 4 == 0 ? 0 : 4 - row_bytes % 4;
	static uint8_t buff[4] = { 0, 0, 0, 0 };
	if (!bitmap || !info || !stream.is_open())
		return false;
	else if (info->img_channels == save_img_channels)
	{
#define SAVE_8_BIT(TYPE) RGBa_t palette[256];															\
		get_most_important_color(palette, 256);															\
		SET_IMG_ARR(palette, 256 * sizeof(RGBa_t));														\
		uint8_t* save_bmp = bitmap + (info->img_height - 1) * (info->img_width * info->img_channels);	\
		for (int i = 0; i < info->img_height; ++i)														\
		{																								\
			TYPE* p = reinterpret_cast<TYPE*>(save_bmp);												\
			for (int j = 0; j < info->img_width; ++j)													\
			{																							\
				uint8_t min_index = 0;																	\
				double min_len = 255.0 * 255.0 * 255.0;													\
				for (int k = 0; k < 256; ++k)															\
				{																						\
					uint8_t r = p[j].red - palette[k].red;												\
					uint8_t g = p[j].green - palette[k].green;											\
					uint8_t b = p[j].blue - palette[k].blue;											\
					double current_len = sqrt(r * r + g * g + b * b);									\
					if (current_len < min_len)															\
					{																					\
						min_index = k;																	\
						min_len = current_len;															\
					}																					\
				}																						\
				SET_IMG_BYTES(min_index);																\
			}																							\
			SET_IMG_ARR(buff, need);																	\
			save_bmp -= info->img_width * info->img_channels;											\
		}
		SAVE_8_BIT(RGB_t)
			return true;
	}
	else if (info->img_channels > save_img_channels)
	{
		SAVE_8_BIT(RGBa_t)
			return true;
	}
	else
		return false;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_24bit_t, info)
{
	LOAD_INFO_BITS
	return info->bit_count == 24 && info->compression == 0;
}

IMPLEMENT_LOAD_FUNCTION(O3DS::img_formats::bmp_24bit_t, bitmap)
{
	if (!info || info->bit_count != 24) return false;
	info->img_channels = sizeof(RGB_t);
	c = info->img_channels;
	RGB_t* arr = new RGB_t[info->img_width * info->img_height];
	const int need = info->img_width * sizeof(RGB_t) % 4 == 0 ? 0 : 4 - (info->img_width * sizeof(RGB_t) % 4);
	const int bias = info->img_width * sizeof(RGB_t) + need;
	for (int i = 0; i < info->img_height; ++i)
	{
		uint8_t buff[4];
		stream.seekg(header.offset_bits + i * bias, std::ios_base::beg);
		GET_IMG_ARR(arr + (info->img_height - 1 - i) * info->img_width, info->img_width * sizeof(RGB_t));
		GET_IMG_ARR(buff, need);
	}
	bitmap = reinterpret_cast<uint8_t*>(arr);
	return bitmap;
}
