#ifndef _OPEN_3D_SCAN_STREAMS
#define _OPEN_3D_SCAN_STREAMS

#include <fstream>
#include <string>

namespace O3DS
{
	namespace streams
	{
		struct istream_interface_t
		{
			virtual istream_interface_t& read(char* bytes, size_t count) noexcept = 0;
			virtual istream_interface_t& seekg(size_t, int mode) noexcept = 0;
			virtual bool is_open() const noexcept = 0;
		};

		struct ostream_interface_t
		{
			virtual ostream_interface_t& write(const char* bytes, size_t count) noexcept = 0;
			virtual ostream_interface_t& seekp(size_t count, int mode) noexcept = 0;
			virtual bool is_open() const noexcept = 0;
		};

		struct std_binaryifstream_t : istream_interface_t
		{
			std::ifstream stream;

			std_binaryifstream_t();
			std_binaryifstream_t(std::wstring path);
			void open(std::wstring path, std::ios_base::openmode mode) noexcept;
			void close() noexcept;
			istream_interface_t& read(char* bytes, size_t count) noexcept override;
			istream_interface_t& seekg(size_t count, int mode) noexcept override;
			bool is_open() const noexcept override;
		};

		struct std_binaryofstream_t : ostream_interface_t
		{
			std::ofstream stream;

			std_binaryofstream_t();
			std_binaryofstream_t(std::wstring path);
			void open(std::wstring path, std::ios_base::openmode mode) noexcept;
			void close() noexcept;
			ostream_interface_t& write(const char* bytes, size_t count) noexcept override;
			ostream_interface_t& seekp(size_t count, int mode) noexcept override;
			bool is_open() const noexcept override;
		};

		struct std_textifstream_t : istream_interface_t
		{
			std::ifstream stream;

			std_textifstream_t();
			std_textifstream_t(std::wstring path);
			void open(std::wstring path, std::ios_base::openmode mode) noexcept;
			void close() noexcept;
			istream_interface_t& read(char* bytes, size_t count) noexcept override;
			istream_interface_t& seekg(size_t count, int mode) noexcept override;
			bool is_open() const noexcept override;
		};

		struct std_textofstream_t : ostream_interface_t
		{
			std::ofstream stream;

			std_textofstream_t();
			std_textofstream_t(std::wstring path);
			void open(std::wstring path, std::ios_base::openmode mode) noexcept;
			void close() noexcept;
			ostream_interface_t& write(const char* bytes, size_t count) noexcept override;
			ostream_interface_t& seekp(size_t count, int mode) noexcept override;
			bool is_open() const noexcept override;
		};

		std::wstring get_extension(const std::wstring& path) noexcept;
		std::wstring get_filename(const std::wstring& path, bool with_extension = true) noexcept;
		std::wstring get_root(const std::wstring& path) noexcept;
		bool path_is_directory(const std::wstring& path) noexcept;
		bool path_is_file(const std::wstring& path) noexcept;
	}
}

#endif