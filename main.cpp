#include <iostream>
#include <cmath>
#include <array>

constexpr double C0 = 299792458; /// Speed of light [metres per second]

enum class side_type
{
  left, bottom, right, top
};

constexpr unsigned int side_to_id (side_type side)
{
  switch (side)
    {
      case side_type::left:   return 0;
      case side_type::bottom: return 1;
      case side_type::right:  return 2;
      case side_type::top:    return 3;
      default:                return 0;
    }

  return 0;
}

constexpr unsigned int get_cell_x (unsigned int cell_id, unsigned int nx)
{
  return cell_id % nx;
}

constexpr unsigned int get_cell_y (unsigned int cell_id, unsigned int nx)
{
  return cell_id / nx;
}

constexpr unsigned int get_cell_id (unsigned int x, unsigned int y, unsigned int nx) { return y * nx + x; }

constexpr bool is_edge_left (unsigned int edge_id) { return edge_id == side_to_id (side_type::left); }
constexpr bool is_edge_bottom (unsigned int edge_id) { return edge_id == side_to_id (side_type::bottom); }
constexpr bool is_edge_right (unsigned int edge_id) { return edge_id == side_to_id (side_type::right); }
constexpr bool is_edge_top (unsigned int edge_id) { return edge_id == side_to_id (side_type::top); }

constexpr bool is_cell_on_left_boundary (unsigned int x) { return x == 0; }
constexpr bool is_cell_on_right_boundary (unsigned int x, unsigned int nx) { return x == nx - 1; }
constexpr bool is_cell_on_bottom_boundary (unsigned int y) { return y == 0; }
constexpr bool is_cell_on_top_boundary (unsigned int y, unsigned int ny) { return y == ny - 1; }

constexpr unsigned int get_neighbor_id (unsigned int cell_id, unsigned int edge_id, unsigned int nx, unsigned int ny)
{
  unsigned int x = get_cell_x (cell_id, nx);
  unsigned int y = get_cell_y (cell_id, nx);

  if (is_edge_left (edge_id))
    return is_cell_on_left_boundary (x) ? cell_id : get_cell_id (x - 1, y, nx);
  else if (is_edge_bottom (edge_id))
    return is_cell_on_bottom_boundary (y) ? cell_id : get_cell_id (x, y - 1, nx);
  else if (is_edge_right (edge_id))
    return is_cell_on_right_boundary (x, nx) ? cell_id : get_cell_id (x + 1, y, nx);
  else if (is_edge_top (edge_id))
    return is_cell_on_top_boundary (y, ny) ? cell_id : get_cell_id (x, y + 1, nx);

  return -1;
}

template<unsigned int nx, unsigned int ny>
constexpr float update_curl_ex (unsigned int cell_id, float dy, const std::array<float, nx * ny> &ez)
{
  const unsigned int neighbor_id = get_neighbor_id (cell_id, side_to_id (side_type::top), nx, ny);
  return (ez[neighbor_id] - ez[cell_id]) / dy;
}

template<unsigned int nx, unsigned int ny>
constexpr float update_curl_ey (unsigned int cell_id, float dx, const std::array<float, nx * ny> &ez)
{
  const unsigned int neighbor_id = get_neighbor_id (cell_id, side_to_id (side_type::right), nx, ny);
  return -(ez[neighbor_id] - ez[cell_id]) / dx;
}

template <unsigned int nx, unsigned int ny>
constexpr float calculate_dt (float dx, float dy)
{
  constexpr float cfl = 0.6;

  return cfl * std::min (dx, dy) / C0;
}

template <unsigned int nx, unsigned int ny>
constexpr float update_curl_h (
  unsigned int cell_id,
  float dx, float dy,
  const std::array<float, nx * ny> hx,
  const std::array<float, nx * ny> hy)
{
  // TODO For now assume that only periodic boundary conditions exist
  const unsigned int left_neighbor_id = get_neighbor_id (cell_id, side_to_id (side_type::left), nx, ny);
  const unsigned int bottom_neighbor_id = get_neighbor_id (cell_id, side_to_id (side_type::bottom), nx, ny);

  return (hy[cell_id] - hy[left_neighbor_id])   / dx
       - (hx[cell_id] - hx[bottom_neighbor_id]) / dy;
}

constexpr float gaussian_pulse (float t, float t_0, float tau)
{
  return std::exp (-(((t - t_0) / tau) * (t - t_0) / tau));
}

constexpr float harmonic_source (float t, float frequency)
{
  const float value = t * frequency;

  if (value != 0)
    return std::cos (value);
  return 0.0;
}

constexpr float step_source (float t)
{
  if (t < 1E-10)
    return 1.0f;
  return 0.0f;
}

constexpr float calculate_source (float t, float frequency)
{
  #if 0
  const float tau = 0.5f / frequency;
  const float t_0 = 6.0f * tau;
  return gaussian_pulse (t, t_0, tau);
  #endif

  return harmonic_source (t, frequency);
  // return step_source (t);
}

template <unsigned int nx, unsigned int ny>
constexpr void fdtd (
  float &t,
  float dt,
  float dx,
  float dy,
  unsigned int steps,
  std::array<float, nx * ny> &er,
  const std::array<float, nx * ny> &mh,

  std::array<float, nx * ny> &hx,
  std::array<float, nx * ny> &hy,
  std::array<float, nx * ny> &ez,
  std::array<float, nx * ny> &dz
  )
{
  constexpr unsigned int n_cells = nx * ny;

  const unsigned int source_cell = get_cell_id(nx / 2, ny / 2, nx);
  const float source_frequency = 1E+9;

  for (unsigned int step = 0; step < steps; step++)
    {
      /// Update h
      for (unsigned int i = 0; i < n_cells; i++)
        {
          hx[i] -= mh[i] * update_curl_ex<nx, ny> (i, dy, ez);
          hy[i] -= mh[i] * update_curl_ey<nx, ny> (i, dx, ez);
        }

      /// Update e
      for (unsigned int i = 0; i < n_cells; i++)
        {
          dz[i] += C0 * dt * update_curl_h<nx, ny> (i, dx, dy, hx, hy);
          dz[source_cell] += calculate_source (t, source_frequency); /// Update source
          ez[i] = dz[i] / er[i];
        }

      t += dt;
    }
}

template <unsigned int nx, unsigned int ny, unsigned int report_steps>
constexpr auto collect_time_steps (unsigned int report_each)
{
  constexpr unsigned int n_cells = nx * ny;
  constexpr float dx = 3.0f / nx;
  constexpr float dy = 3.0f / ny;
  constexpr double dt = calculate_dt<nx, ny> (dx, dy);
  std::array<float, n_cells> er {};
  std::array<float, n_cells> hr {};
  std::array<float, n_cells> mh {};

  std::array<float, n_cells> hx {};
  std::array<float, n_cells> hy {};
  std::array<float, n_cells> ez {};
  std::array<float, n_cells> dz {};

  for (unsigned int i = 0; i < n_cells; i++)
    {
      er[i] = 1.0;
      hr[i] = 1.0;

      hx[i] = 0.0;
      hy[i] = 0.0;
      ez[i] = 0.0;
      dz[i] = 0.0;
    }

  for (unsigned int i = 0; i < n_cells; i++)
    mh[i] = C0 * dt / hr[i];

  float t = 0.0;

  std::array<std::array<float, nx * ny>, report_steps> report {};
  for (unsigned int report_step = 0; report_step < report_steps; report_step++)
    {
      fdtd<nx, ny> (t, dt, dx, dy, report_each, er, mh, hx, hy, ez, dz);

      for (unsigned int i = 0; i < n_cells; i++)
        report[report_step][i] = ez[i];
    }
  return report;
}

template <class float_type>
void write_vtk (
  const std::string &filename,
  const float_type dx,
  const float_type dy,
  const unsigned int nx,
  const unsigned int ny,
  const float_type *e)
{
  FILE * f = fopen (filename.c_str (), "w");

  fprintf (f, "# vtk DataFile Version 3.0\n");
  fprintf (f, "vtk output\n");
  fprintf (f, "ASCII\n");
  fprintf (f, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (f, "POINTS %u double\n", nx * ny * 4);

  for (unsigned int j = 0; j < ny; j++)
    {
      for (unsigned int i = 0; i < nx; i++)
        {
          fprintf (f, "%lf %lf 0.0\n", dx * (i + 0), dy * (j + 0) );
          fprintf (f, "%lf %lf 0.0\n", dx * (i + 1), dy * (j + 0) );
          fprintf (f, "%lf %lf 0.0\n", dx * (i + 1), dy * (j + 1) );
          fprintf (f, "%lf %lf 0.0\n", dx * (i + 0), dy * (j + 1) );
        }
    }

  fprintf (f, "CELLS %u %u\n", nx * ny, nx * ny * 5);

  for (unsigned int j = 0; j < ny; j++)
    {
      for (unsigned int i = 0; i < nx; i++)
        {
          const unsigned int point_offset = (j * nx + i) * 4;
          fprintf (f, "4 %u %u %u %u\n", point_offset + 0, point_offset + 1, point_offset + 2, point_offset + 3);
        }
    }

  fprintf (f, "CELL_TYPES %u\n", nx * ny);
  for (unsigned int i = 0; i < nx * ny; i++)
    fprintf (f, "9\n");

  fprintf (f, "CELL_DATA %u\n", nx * ny);
  fprintf (f, "SCALARS Ez double 1\n");
  fprintf (f, "LOOKUP_TABLE default\n");

  for (unsigned int i = 0; i < nx * ny; i++)
    fprintf (f, "%lf\n", e[i]);

  fclose (f);
}

int main ()
{
  constexpr unsigned int nx = 45;
  constexpr unsigned int ny = 45;
  constexpr auto reports = collect_time_steps<nx, ny, 8> (10);

  for (unsigned int report = 0; report < reports.size (); report++)
    write_vtk<float> ("output_" + std::to_string (report) + ".vtk", 3.0 / nx, 3.0 / ny, nx, ny, reports[report].data ());

  return 0;
}
