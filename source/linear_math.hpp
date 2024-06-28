#ifndef _OPEN_3D_SCAN_VECTOR_MATRIX
#define _OPEN_3D_SCAN_VECTOR_MATRIX

#include <iostream>
#include <array>
#include <algorithm>

namespace O3DS
{
	namespace math
	{
		template <typename type_t, size_t dim>
		class vector final
		{
		public:
			std::array<type_t, dim> _data;

			template <typename elem_t, typename ...args_t>
			void _fill(size_t index, const elem_t& first, const args_t&... other) noexcept
			{
				static_assert(sizeof...(other) < dim, "list overflow");
				_data[index] = static_cast<type_t>(first);
				_fill(index + 1, other...);
			}

			template <typename elem_t>
			void _fill(size_t index, const elem_t& last) noexcept { _data[index] = static_cast<type_t>(last); }

			size_t _hash(const std::string& key) const noexcept
			{
				if (key == "x") return 0;
				if (key == "y") return 1;
				if (key == "z") return 2;
				if (key == "w") return 3;
				return _data.size();
			}

			friend class vector;
		public:
			using value_type = type_t;

			using iterator = decltype(_data.begin());
			using const_iterator = decltype(_data.cbegin());
			using reverse_iterator = decltype(_data.rbegin());
			using const_reverse_iterator = decltype(_data.crbegin());

			vector() noexcept : _data() {}

			explicit vector(const type_t& value) noexcept
			{
				for (auto i = _data.begin(); i != _data.end(); *i = value, ++i);
			}

			template <typename ..._args>
			vector(const _args&... values) noexcept : _data() { _fill(0, values...); }

			template <typename other_type_t>
			vector(const vector<other_type_t, dim>& other) noexcept
			{
				std::copy(other.cbegin(), other.cend(), _data.begin());
			}

			template <typename other_type_t>
			vector<type_t, dim>& operator=(const vector<other_type_t, dim>& other) noexcept
			{
				std::copy(other.cbegin(), other.cend(), _data.begin());
				return *this;
			}

			vector(const vector<type_t, dim>&) = default;
			vector<type_t, dim>& operator=(const vector<type_t, dim>&) = default;

			vector(vector<type_t, dim>&&) = default;
			vector& operator=(vector<type_t, dim>&&) = default;

			iterator				begin()		noexcept { return _data.begin(); }
			iterator				end()		noexcept { return _data.end(); }
			const_iterator			cbegin()	const noexcept { return _data.cbegin(); }
			const_iterator			cend()		const noexcept { return _data.cend(); }
			reverse_iterator		rbegin()	noexcept { return _data.rbegin(); }
			reverse_iterator		rend()		noexcept { return _data.rend(); }
			const_reverse_iterator	crbegin()	const noexcept { return _data.crbegin(); }
			const_reverse_iterator	crend()		const noexcept { return _data.crend(); }

			static constexpr size_t dimension = dim;
			const type_t* ptr() const noexcept { return _data.data(); }

			type_t& operator[](size_t index) { return _data[index]; }
			const type_t& operator[](size_t index) const { return _data[index]; }

			type_t& operator[](const std::string& key) { return _data[_hash(key)]; }
			const type_t& operator[](const std::string& key) const { return _data[_hash(key)]; }

			vector<type_t, dim>& operator+=(const vector<type_t, dim>& other) noexcept
			{
				for (std::size_t i = 0; i < dim; ++i)
					_data[i] += other._data[i];
				return *this;
			}

			vector<type_t, dim>& operator-=(const vector<type_t, dim>& other) noexcept
			{
				for (std::size_t i = 0; i < dim; ++i)
					_data[i] -= other._data[i];
				return *this;
			}
		};

		using vec2 = vector<double, 2>;
		using vec3 = vector<double, 3>;
		using vec4 = vector<double, 4>;

		using vec2f = vector<float, 2>;
		using vec3f = vector<float, 3>;
		using vec4f = vector<float, 4>;

		using vec2i = vector<int, 2>;
		using vec3i = vector<int, 3>;
		using vec4i = vector<int, 4>;

		using size2 = vector<size_t, 2>;
		using size3 = vector<size_t, 3>;

		template <typename type_t, size_t rows, size_t columns>
		class matrix final
		{
		private:
			std::array<type_t, rows* columns> _data;

			template <typename elem_t, typename ...args_t>
			void _fill(size_t index, const elem_t& first, const args_t&... other) noexcept
			{
				static_assert(sizeof...(other) < rows * columns, "list overflow");
				_data[index] = static_cast<type_t>(first);
				_fill(index + 1, other...);
			}

			template <typename elem_t>
			void _fill(size_t index, const elem_t& last) noexcept { _data[index] = static_cast<type_t>(last); }

			size_t _hash(const std::string& key) const noexcept
			{
				if (key == "x")  return 0;
				if (key == "xy") return 1;
				if (key == "xz") return 2;
				if (key == "xw") return 3;
				if (key == "yx") return 0 + columns;
				if (key == "y")  return 1 + columns;
				if (key == "yz") return 2 + columns;
				if (key == "yw") return 3 + columns;
				if (key == "zx") return 0 + 2 * columns;
				if (key == "zy") return 1 + 2 * columns;
				if (key == "z")  return 2 + 2 * columns;
				if (key == "zw") return 3 + 2 * columns;
				if (key == "wx") return 0 + 3 * columns;
				if (key == "wy") return 1 + 3 * columns;
				if (key == "wz") return 2 + 3 * columns;
				if (key == "w")  return 3 + 3 * columns;
				return _data.size();
			}

			friend class matrix;
		public:
			using value_type = type_t;

			using iterator = decltype(_data.begin());
			using const_iterator = decltype(_data.cbegin());
			using reverse_iterator = decltype(_data.rbegin());
			using const_reverse_iterator = decltype(_data.crbegin());

			matrix() noexcept : _data()
			{
				if constexpr (rows == columns)
					for (size_t i = 0; i < rows; _data[i + i * columns] = static_cast<type_t>(1.0), ++i);
			}

			explicit matrix(const type_t& value) noexcept : _data()
			{
				if constexpr (rows == columns)
					for (size_t i = 0; i < rows; _data[i + i * columns] = value, ++i);
				else
					for (size_t i = 0; i < rows * columns; _data[i] = value, ++i);
			}

			explicit matrix(const vector<type_t, rows>& values) noexcept : _data()
			{
				static_assert(rows == columns, "matrix is not square");
				for (size_t i = 0; i < rows; ++i)
					_data[i + i * columns] = values[i];
			}

			template <typename ...args_t>
			matrix(const args_t&... values) noexcept : _data() { _fill(0, values...); }

			template <typename other_type_t>
			matrix(const matrix<other_type_t, rows, columns>& other) noexcept
			{
				std::copy(other.cbegin(), other.cend(), _data.begin());
			}

			template <typename other_type_t>
			matrix<type_t, rows, columns>& operator=(const matrix<other_type_t, rows, columns>& other) noexcept
			{
				std::copy(other.cbegin(), other.cend(), _data.begin());
				return *this;
			}

			matrix(const matrix<type_t, rows, columns>&) = default;
			matrix<type_t, rows, columns>& operator=(const matrix<type_t, rows, columns>&) = default;

			matrix(matrix<type_t, rows, columns>&&) = default;
			matrix& operator=(matrix<type_t, rows, columns>&&) = default;

			iterator				begin()		noexcept { return _data.begin(); }
			iterator				end()		noexcept { return _data.end(); }
			const_iterator			cbegin()	const noexcept { return _data.cbegin(); }
			const_iterator			cend()		const noexcept { return _data.cend(); }
			reverse_iterator		rbegin()	noexcept { return _data.rbegin(); }
			reverse_iterator		rend()		noexcept { return _data.rend(); }
			const_reverse_iterator	crbegin()	const noexcept { return _data.crbegin(); }
			const_reverse_iterator	crend()		const noexcept { return _data.crend(); }

			static constexpr size_t row_count = rows;
			static constexpr size_t column_count = columns;
			const type_t* ptr() const noexcept { return _data.data(); }

			type_t* operator[](size_t index) { return _data.data() + index * columns; }
			const type_t* operator[](size_t index) const { return _data.data() + index * columns; }

			type_t& operator[](const std::string& key) { return _data[_hash(key)]; }
			const type_t& operator[](const std::string& key) const { return _data[_hash(key)]; }

			matrix<type_t, rows, columns>& operator+=(const matrix<type_t, rows, columns>& other) noexcept
			{
				for (size_t i = 0; i < rows * columns; ++i)
					_data[i] += other._data[i];
				return *this;
			}

			matrix<type_t, rows, columns>& operator-=(const matrix<type_t, rows, columns>& other) noexcept
			{
				for (size_t i = 0; i < rows * columns; ++i)
					_data[i] -= other._data[i];
				return *this;
			}
		};

		using mat2 = matrix<double, 2, 2>;
		using mat3 = matrix<double, 3, 3>;
		using mat4 = matrix<double, 4, 4>;

		using mat2f = matrix<float, 2, 2>;
		using mat3f = matrix<float, 3, 3>;
		using mat4f = matrix<float, 4, 4>;

		using mat2i = matrix<int, 2, 2>;
		using mat3i = matrix<int, 3, 3>;
		using mat4i = matrix<int, 4, 4>;

		template <typename type_t, size_t dim>
		static vector<type_t, dim> operator+(const vector<type_t, dim>& a, const vector<type_t, dim>& b) noexcept
		{
			vector<type_t, dim> result;
			for (size_t i = 0; i < dim; ++i)
				result[i] = a[i] + b[i];
			return result;
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim> operator-(const vector<type_t, dim>& a, const vector<type_t, dim>& b) noexcept
		{
			vector<type_t, dim> result;
			for (size_t i = 0; i < dim; ++i)
				result[i] = a[i] - b[i];
			return result;
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim> operator*(const vector<type_t, dim>& a, const type_t& scalar) noexcept
		{
			vector<type_t, dim> result;
			for (size_t i = 0; i < dim; ++i)
				result[i] = a[i] * scalar;
			return result;
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim> operator*(const type_t& scalar, const vector<type_t, dim>& a) noexcept
		{
			vector<type_t, dim> result;
			for (size_t i = 0; i < dim; ++i)
				result[i] = a[i] * scalar;
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static matrix<type_t, rows, columns> operator+(
			const matrix<type_t, rows, columns>& A, const matrix<type_t, rows, columns>& B) noexcept
		{
			matrix<type_t, rows, columns> result;
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					result[i][j] = A[i][j] + B[i][j];
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static matrix<type_t, rows, columns> operator-(
			const matrix<type_t, rows, columns>& A, const matrix<type_t, rows, columns>& B) noexcept
		{
			matrix<type_t, rows, columns> result;
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					result[i][j] = A[i][j] - B[i][j];
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns, size_t other_columns>
		static matrix<type_t, rows, columns> operator*(
			const matrix<type_t, rows, columns>& A, const matrix<type_t, columns, other_columns>& B) noexcept
		{
			matrix<type_t, rows, other_columns> result(0.0);
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < other_columns; ++j)
				{
					for (size_t k = 0; k < columns; ++k)
						result[i][j] += A[i][k] * B[k][j];
				}
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static vector<type_t, columns> operator*(
			const matrix<type_t, rows, columns>& A, const vector<type_t, columns>& v) noexcept
		{
			vector<type_t, columns> result(0.0);
			for (size_t j = 0; j < columns; ++j)
			{
				for (size_t i = 0; i < rows; ++i)
					result[i] += v[j] * A[i][j];
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static matrix<type_t, rows, columns> operator*(
			const matrix<type_t, rows, columns>& A, const type_t& scalar) noexcept
		{
			matrix<type_t, rows, columns> result;
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					result[i][j] = A[i][j] * scalar;
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static matrix<type_t, rows, columns> operator*(
			const type_t& scalar, const matrix<type_t, rows, columns>& A) noexcept
		{
			matrix<type_t, rows, columns> result;
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
					result[i][j] = A[i][j] * scalar;
			}
			return result;
		}

		template <typename type_t, size_t dim>
		static bool operator==(const vector<type_t, dim>& a, const vector<type_t, dim>& b) noexcept
		{
			for (size_t i = 0; i < dim; ++i)
			{
				if (a[i] != b[i])
					return false;
			}
			return true;
		}

		template <typename type_t, size_t dim>
		static bool operator!=(const vector<type_t, dim>& a, const vector<type_t, dim>& b) noexcept { return !(a == b); }

		template <typename type_t, size_t rows, size_t columns>
		static bool operator==(const matrix<type_t, rows, columns>& A,
			const matrix<type_t, rows, columns>& B) noexcept
		{
			for (size_t i = 0; i < rows; ++i)
			{
				for (size_t j = 0; j < columns; ++j)
				{
					if (A[i][j] != B[i][j])
						return false;
				}
			}
			return true;
		}

		template <typename type_t, size_t rows, size_t columns>
		static bool operator!=(const matrix<type_t, rows, columns>& A,
			const matrix<type_t, rows, columns>& B) noexcept { return !(A == B); }

		template <typename type_t, size_t dim>
		static double sum(const vector<type_t, dim>& v) noexcept
		{
			double result = 0.0;
			for (int i = 0; i < dim; ++i)
				result += v[i];
			return result;
		}

		template <typename type_t, size_t dim>
		static double length(const vector<type_t, dim>& vec) noexcept
		{
			double result = 0.0;
			for (size_t i = 0; i < dim; result += vec[i] * vec[i], ++i);
			return sqrt(result);
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim> normalize(const vector<type_t, dim>& vec)
		{
			vector<type_t, dim> result;
			double len = length(vec);
			for (size_t i = 0; i < dim; result[i] = vec[i] / len, ++i);
			return result;
		}

		template <typename type_t, size_t current_rows, size_t current_columns>
		static double det(const matrix<type_t, current_rows, current_columns>& mat) noexcept
		{
			if constexpr (current_rows != current_columns)
				return 0.0;
			constexpr size_t n = current_rows;
			if constexpr (n <= 1)
				return mat[0][0];
			else if constexpr (n == 2)
				return (mat["x"] * mat["y"]) - (mat["xy"] * mat["yx"]);
			else
			{
				double result = 0.0;
				for (size_t k = 0; k < n; k++)
				{
					matrix<type_t, current_rows - 1, current_columns - 1> minor;
					for (size_t i = 0; i < n - 1; i++)
					{
						bool column = false;
						for (size_t j = 0; j < n - 1; j++)
						{
							if (j == k) column = true;
							minor[i][j] = column ? mat[i + 1][j + 1] : mat[i + 1][j];
						}
					}
					result += pow(-1.0, k) * det(minor) * mat[0][k];
				}
				return result;
			}
		}

		template <typename type_t, size_t dim>
		static double dot(const vector<type_t, dim>& a, const vector<type_t, dim>& b) noexcept
		{
			double result = 0;
			for (size_t i = 0; i < dim; ++i)
				result += a[i] * b[i];
			return result;
		}

		template <typename type_t>
		static vector<type_t, 3> cross(const vector<type_t, 3>& a, const vector<type_t, 3>& b) noexcept
		{
			matrix<type_t, 2, 2> mat_i = {
				a["y"], b["y"],
				a["z"], b["z"]
			};
			matrix<type_t, 2, 2> mat_j = {
				a["x"], b["x"],
				a["z"], b["z"]
			};
			matrix<type_t, 2, 2> mat_k = {
				a["x"], b["x"],
				a["y"], b["y"]
			};
			return { det(mat_i), -det(mat_j), det(mat_k) };
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim - 1> low(const vector<type_t, dim>& v) noexcept
		{
			vector<type_t, dim - 1> result;
			for (size_t i = 0; i < dim - 1; ++i)
				result[i] = v[i];
			return result;
		}

		template <typename type_t, size_t dim>
		static matrix<type_t, dim - 1, dim - 1> low(const matrix<type_t, dim, dim>& M) noexcept
		{
			matrix<type_t, dim - 1, dim - 1> result;
			for (size_t i = 0; i < dim - 1; ++i)
			{
				for (size_t j = 0; j < dim - 1; ++j)
					result[i][j] = M[i][j];
			}
			return result;
		}

		template <typename type_t, size_t dim>
		static vector<type_t, dim + 1> high(const vector<type_t, dim>& v, const type_t& lost) noexcept
		{
			vector<type_t, dim + 1> result;
			for (size_t i = 0; i < dim; ++i)
				result[i] = v[i];
			result[dim] = lost;
			return result;
		}

		template <typename type_t, size_t dim>
		static matrix<type_t, dim + 1, dim + 1> high(
			const matrix<type_t, dim, dim>& M, const type_t& lost) noexcept
		{
			matrix<type_t, dim + 1, dim + 1> result(0.0);
			for (size_t i = 0; i < dim; ++i)
			{
				for (size_t j = 0; j < dim; ++j)
					result[i][j] = M[i][j];
			}
			result[dim][dim] = lost;
			return result;
		}

		template <typename type_t, size_t dim>
		static matrix<type_t, dim, dim> transpose(const matrix<type_t, dim, dim>& M) noexcept
		{
			matrix<type_t, dim, dim> result;
			for (size_t i = 0; i < dim; ++i)
			{
				for (size_t j = 0; j < dim; ++j)
					result[i][j] = M[j][i];
			}
			return result;
		}

		template <typename type_t, size_t rows, size_t columns>
		static matrix<type_t, rows, columns> change(
			const matrix<type_t, rows, columns>& M, const matrix<type_t, rows, columns>& basis) noexcept
		{
			return tranpose(basis) * M * basis;
		}

		static mat3 dcm_psi(double psi) noexcept
		{
			return {
				cos(psi), 0.0, sin(psi),
				0.0, 1.0, 0.0,
				-sin(psi), 0.0, cos(psi)
			};
		}

		static mat3 dcm_tetta(double tetta) noexcept
		{
			return {
				1.0, 0.0, 0.0,
				0.0, sin(tetta), cos(tetta),
				0.0, cos(tetta), -sin(tetta)
			};
		}

		static mat3 dcm_gamma(double gamma) noexcept
		{
			return {
				cos(gamma), sin(gamma), 0.0,
				sin(gamma), cos(gamma), 0.0,
				0.0, 0.0, 1.0
			};
		}

		// psi, tetta, gamma
		static mat3 dcm_g21(vec3 euler) noexcept
		{
			return dcm_gamma(euler["z"]) * dcm_tetta(euler["y"]) * dcm_psi(euler["x"]);
		}

		// psi, tetta, gamma
		static mat3 dcm_12g(vec3 euler) noexcept
		{
			return transpose(dcm_g21(euler));
		}
	}
}

#endif
