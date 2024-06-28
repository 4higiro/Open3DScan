#ifndef _OPEN_3D_SCAN_NUMERICS
#define _OPEN_3D_SCAN_NUMERICS

#include "linear_math.hpp"
#include <math.h>

namespace O3DS
{
	namespace math
	{
		static constexpr double pi = 3.14159265359;
		static constexpr size_t inf = 0b11111111'11111111'11111111'11111111'11111111'11111111'11111111'11111111;
		static constexpr long long positive_inf = 0b11111111'11111111'11111111'11111111'11111111'11111111'11111111'1111111;
		static constexpr long long negative_inf = -0b11111111'11111111'11111111'11111111'11111111'11111111'11111111'1111111;

		static constexpr double rad(double deg_angle) noexcept { return deg_angle * (pi / 180.0); }
		static constexpr double deg(double rad_angle) noexcept { return rad_angle * 180.0 / pi; }

		static constexpr int round(long double number) noexcept
		{
			if (number > 0.0)
				return number - static_cast<int>(number) > 0.5 ? static_cast<int>(number + 1) : static_cast<int>(number);
			else
				return number - static_cast<int>(number) < -0.5 ? static_cast<int>(number - 1) : static_cast<int>(number);
		}

		static constexpr double sqrt_2pi = 2.50662827463;

		static double gauss(double t, double sigma) noexcept
		{
			return (1.0 / (sigma * sqrt_2pi)) * exp(-t * t / (2.0 * sigma * sigma));
		}

		static double d_dt_gauss(double t, double sigma) noexcept
		{
			return (-t / (sigma * sigma * sigma * sqrt_2pi)) * exp(-t * t / (2.0 * sigma * sigma));
		}

		static double d2_dt2_gauss(double t, double sigma) noexcept
		{
			return ((t * t - sigma * sigma) * exp(-t * t / (2.0 * sigma * sigma))) /
				(sigma * sigma * sigma * sigma * sigma * sqrt_2pi);
		}

		template <size_t dim>
		static vector<int, dim> round(const vector<float, dim>& v) noexcept
		{
			vector<int, dim> result;
			for (int i = 0; i < dim; ++i)
				result[i] = round(v[i]);
			return result;
		}

		template <size_t dim>
		static vector<int, dim> round(const vector<double, dim>& v) noexcept
		{
			vector<int, dim> result;
			for (int i = 0; i < dim; ++i)
				result[i] = round(v[i]);
			return result;
		}

		template <size_t dim>
		static vector<int, dim> round(const vector<long double, dim>& v) noexcept
		{
			vector<int, dim> result;
			for (int i = 0; i < dim; ++i)
				result[i] = round(v[i]);
			return result;
		}

		template <typename type_t, size_t dim> // kramer method
		static vector<type_t, dim> get_slau_solution(const matrix<type_t, dim, dim>& koeffs, const vector<type_t, dim>& free_terms) noexcept
		{
			vector<type_t, dim> result;
			double delta = det(koeffs);
			if (!delta) return result;
			for (int k = 0; k < dim; k++)
			{
				matrix<type_t, dim, dim> delta_matrix;
				for (int i = 0; i < dim; i++)
				{
					for (int j = 0; j < dim; j++)
					{
						if (j == k)
							delta_matrix[i][j] = free_terms[i];
						else
							delta_matrix[i][j] = koeffs[i][j];
					}
				}
				result[k] = det(delta_matrix) / delta;
			}
			return result;
		}

		template <typename type_t, size_t dim>
		static matrix<type_t, dim, dim> inverse(const matrix<type_t, dim, dim>& M) noexcept
		{
			matrix<type_t, dim, dim> result;
			for (int i = 0; i < dim; ++i)
			{
				vector<type_t, dim> free_terms;
				free_terms[i] = 1;
				vector<type_t, dim> v = get_slau_solution(M, free_terms);
				for (int j = 0; j < dim; ++j)
					result[i][j] = v[j];
			}
			return transpose(result);
		}

		static vec2 get_eigenvalues(const mat2& mat) noexcept
		{
			double D = (mat["x"] + mat["y"]) * (mat["x"] + mat["y"]) - 4.0 * (mat["x"] * mat["y"] - mat["xy"] * mat["yx"]);
			if (D < 0) return { negative_inf, positive_inf };
			double l1 = (mat["x"] + mat["y"] - sqrt(D)) / 2.0;
			double l2 = (mat["x"] + mat["y"] + sqrt(D)) / 2.0;
			return { l1, l2 };
		}

		static vec2 get_eigenvector(const mat2& mat, double lambda, double noise_threshold = 1e-12) noexcept
		{
			if (abs(mat["y"] - lambda) < noise_threshold && abs(mat["x"] - lambda) < noise_threshold)
				return { negative_inf, positive_inf };
			else if (abs(mat["y"] - lambda) < noise_threshold)
			{
				double y = 1.0;
				double x = -mat["xy"] / (mat["x"] - lambda);
				vec2 v = { x, y };
				v = v * (1.0 / sqrt(x * x + y * y)) * (lambda < 0.0 ? -1.0 : 1.0);
				return v;
			}
			else
			{
				double x = 1.0;
				double y = -mat["yx"] / (mat["y"] - lambda);
				vec2 v = { x, y };
				v = v * (1.0 / sqrt(x * x + y * y)) * (lambda < 0.0 ? -1.0 : 1.0);
				return v;
			}
		}

		template <typename function_t>
		static double partial_derivative_x(function_t f, double x, double y, double h = 0.01) noexcept
		{
			return (f(x + h, y) - f(x, y)) / h;
		}

		template <typename function_t>
		static double partial_derivative_y(function_t f, double x, double y, double h = 0.01) noexcept
		{
			return (f(x, y + h) - f(x, y + h)) / h;
		}

		template <typename function_t>
		static math::vec2 gradient(function_t f, double x, double y, double h = 0.01) noexcept
		{
			return { partial_derivative_x(f, x, y, h), partial_derivative_y(f, x, y, h) };
		}
	}
}

#endif
