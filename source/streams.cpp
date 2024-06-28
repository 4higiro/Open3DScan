#include "streams.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

O3DS::streams::std_binaryifstream_t::std_binaryifstream_t() = default;

O3DS::streams::std_binaryifstream_t::std_binaryifstream_t(std::wstring path)
	: stream(path, std::ios_base::binary) {}

void O3DS::streams::std_binaryifstream_t::open(std::wstring path, std::ios_base::openmode mode) noexcept
{
	stream.open(path, mode);
}

void O3DS::streams::std_binaryifstream_t::close() noexcept
{
	stream.close();
}

O3DS::streams::istream_interface_t& O3DS::streams::std_binaryifstream_t::read(char* bytes, size_t count) noexcept
{
	stream.read(bytes, count);
	return *this;
}

O3DS::streams::istream_interface_t& O3DS::streams::std_binaryifstream_t::seekg(size_t count, int mode) noexcept
{
	stream.seekg(count, mode);
	return *this;
}

bool O3DS::streams::std_binaryifstream_t::is_open() const noexcept
{
	return stream.is_open();
}

O3DS::streams::std_binaryofstream_t::std_binaryofstream_t() = default;

O3DS::streams::std_binaryofstream_t::std_binaryofstream_t(std::wstring path)
	: stream(path, std::ios_base::binary) {}

void O3DS::streams::std_binaryofstream_t::open(std::wstring path, std::ios_base::openmode mode) noexcept
{
	stream.open(path, mode);
}

void O3DS::streams::std_binaryofstream_t::close() noexcept
{
	stream.close();
}

O3DS::streams::ostream_interface_t& O3DS::streams::std_binaryofstream_t::write(const char* bytes, size_t count) noexcept
{
	stream.write(bytes, count);
	return *this;
}

O3DS::streams::ostream_interface_t& O3DS::streams::std_binaryofstream_t::seekp(size_t count, int mode) noexcept
{
	stream.seekp(count, mode);
	return *this;
}

bool O3DS::streams::std_binaryofstream_t::is_open() const noexcept
{
	return stream.is_open();
}

O3DS::streams::std_textifstream_t::std_textifstream_t() = default;

O3DS::streams::std_textifstream_t::std_textifstream_t(std::wstring path)
	: stream(path) {}

void O3DS::streams::std_textifstream_t::open(std::wstring path, std::ios_base::openmode mode) noexcept
{
	stream.open(path, mode);
}

void O3DS::streams::std_textifstream_t::close() noexcept
{
	stream.close();
}

O3DS::streams::istream_interface_t& O3DS::streams::std_textifstream_t::read(char* bytes, size_t count) noexcept
{
	std::string text;
	stream >> text;
	for (int i = 0; i < count && i < text.size(); ++i)
		bytes[i] = text[i];
	return *this;
}

O3DS::streams::istream_interface_t& O3DS::streams::std_textifstream_t::seekg(size_t count, int mode) noexcept
{
	stream.seekg(count, mode);
	return *this;
}

bool O3DS::streams::std_textifstream_t::is_open() const noexcept
{
	return stream.is_open();
}

O3DS::streams::std_textofstream_t::std_textofstream_t() = default;

O3DS::streams::std_textofstream_t::std_textofstream_t(std::wstring path)
	: stream(path) {}

void O3DS::streams::std_textofstream_t::open(std::wstring path, std::ios_base::openmode mode) noexcept
{
	stream.open(path, mode);
}

void O3DS::streams::std_textofstream_t::close() noexcept
{
	stream.close();
}

O3DS::streams::ostream_interface_t& O3DS::streams::std_textofstream_t::write(const char* bytes, size_t count) noexcept
{
	stream << bytes;
	return *this;
}

O3DS::streams::ostream_interface_t& O3DS::streams::std_textofstream_t::seekp(size_t count, int mode) noexcept
{
	stream.seekp(count, mode);
	return *this;
}

bool O3DS::streams::std_textofstream_t::is_open() const noexcept
{
	return stream.is_open();
}

std::wstring O3DS::streams::get_extension(const std::wstring& path) noexcept
{
	size_t pos1 = path.find_last_of(L".");
	size_t pos2 = path.find_last_of(L"\\");
	return pos2 < pos1 || pos2 == std::wstring::npos ? path.substr(pos1) : L"";
}

std::wstring O3DS::streams::get_filename(const std::wstring& path, bool with_extension) noexcept
{
	size_t pos1 = path.find_last_of(L"\\");
	if (with_extension)
		return path.substr(pos1 + 1);
	size_t pos2 = path.find_last_of(L".");
	if (pos2 == std::wstring::npos)
		return path.substr(pos1 + 1);
	return pos2 > pos1 ? path.substr(pos1 + 1, pos2 - pos1 - 1) : path.substr(pos1 + 1);
}

std::wstring O3DS::streams::get_root(const std::wstring& path) noexcept
{
	size_t pos = path.find_last_of(L"\\");
	return pos == std::wstring::npos ? L"" : path.substr(0, pos);
}

bool O3DS::streams::path_is_directory(const std::wstring& path) noexcept
{
	return path[path.length() - 1] == L'\\';
}

bool O3DS::streams::path_is_file(const std::wstring& path) noexcept
{
	return !path_is_directory(path);
}
