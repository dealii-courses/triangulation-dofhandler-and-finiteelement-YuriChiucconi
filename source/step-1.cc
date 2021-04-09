/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 * 
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

using namespace dealii;

class MyTuple
{
  std::tuple<int, int, int> Tuple;

public:
  MyTuple(Triangulation<2> &T)
    : Tuple(std::make_tuple(T.n_levels(), T.n_cells(), T.n_active_cells())) {}

  void print()
  {
    std::cout << "Number of level:\t" << std::get<0>(this->Tuple) << "\n"
              << "Number of cells:\t" << std::get<1>(this->Tuple) << "\n"
              << "Number of active cells:\t" << std::get<2>(this->Tuple) << std::endl;
  }
};

MyTuple helper(Triangulation<2> &T) {return MyTuple(T);}

void
circle_grid()
{
  const Point<2> center(0, 0);
  const double   radius = 1.0;

  Triangulation<2> triangulation;
  GridGenerator::hyper_ball(triangulation, center, radius, false);
  triangulation.refine_global(2);

  std::ofstream out("circle_grid.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to circle_grid.svg" << std::endl;

  Triangulation<2> triangulation2;
  GridGenerator::hyper_ball(triangulation2, center, radius, true);
  triangulation2.refine_global(2);

  std::ofstream out2("circle_grid2.svg");
  GridOut       grid_out2;
  grid_out2.write_svg(triangulation2, out2);
  std::cout << "Grid written to circle_grid2.svg" << std::endl;
}

void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);

  std::cout << "Number of original vertices: " << triangulation.n_vertices()
            << std::endl;

  triangulation.refine_global(4);

  std::cout << "Numbe  r of vertices after 4 refinmentss: "
            << triangulation.n_vertices() << std::endl;

  std::ofstream out("grid-1.vtk");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.vtk" << std::endl;
  
  MyTuple t = helper(triangulation);
  t.print();
}



void
second_grid()
{
  Triangulation<2> triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);

  // triangulation.reset_manifold(0);
  /*
  questo comando fa si che la triangolazione passi dal descrivere 
  il Manifold che aveva in precedenza (in questo caso una hyper shell)
  al descrivere un FlatManifold
  */


  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.vtk" << std::endl;
  
  MyTuple t = helper(triangulation);
  t.print();
}

void
third_grid()
{
  Triangulation<2> triangulation;
  GridGenerator::hyper_L(triangulation, 0., 0.9);
  /*triangulation.refine_global(1);*/
  const Point<2> corner(0.45, 0.45);

  for (unsigned int step = 0; step < 6; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          Point<2>     cell_center          = cell->center();
          const double distance_from_corner = corner.distance(cell_center);

          if (distance_from_corner < 1. / 3)
            {
              cell->set_refine_flag();
            }
        }
      triangulation.execute_coarsening_and_refinement();
    }

  std::ofstream out("grid-3.vtk");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-3.vtk" << std::endl;
  
  MyTuple t = helper(triangulation);
  t.print();
}


int
main()
{
  first_grid();
  second_grid();
  third_grid();

  circle_grid();
}