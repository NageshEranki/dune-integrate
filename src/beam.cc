// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

//Grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gnuplot.hh>

//Geometry includes
#include <dune/geometry/quadraturerules.hh>

//Function includes
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

//Numerical stuff
#include <dune/common/fvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>


// My class for the Hermite Cubic shape functions on the real line
//
// WARNING: Use only with one-dimensional domains!
template<class DomainType>
class HermiteFiniteElement
{
	public:
		HermiteFiniteElement() {}

		void evaluateFunction(const DomainType& x, std::vector<Dune::FieldVector<double,1>>& out)
		{
			out.resize(size());
			// out = 0; 		 // Check if reinitialisation is necessary!

			// Is it necessary to index the input x??

			out[0] = (x-1)*(x-1)*(2*x+1);

			out[1] = x*(x-1)*(x-1);

			out[2] = x*x*(3-2*x);

			out[3] = x*x*(x-1);
		}

		void evaluateJacobian()
		{
		
		}

		unsigned int size()
		{
			return 4;
		}

};

int main(int argc, char** argv)
{
  try
  {
	//Start!
	
	// Dimension of domain.
	const int dim = 1;

	typedef Dune::BCRSMatrix<double> MatrixType;
	typedef Dune::BlockVector<double> VectorType;
	
	// Structured grid from YaspGrid
	typedef Dune::YaspGrid<dim> GridType;
	GridType grid({1.0},{10});

	//Obtain view into grid
	typedef GridType::LeafGridView GridView;
	GridView gv = grid.leafGridView();

	//Pointer to first element in the grid
	auto it = gv.begin<0>();

	//Obtain Quadrature rules from geometry type
	const auto& quad = Dune::QuadratureRules<double,dim>::rule(it->geometry().type(),2);

	// We are not using the Lagrange shape functions for the beam problem
	//typedef Dune::Functions::LagrangeBasis<GridView,2> BasisType;
	//BasisType basis(gv);
	
	// Create local views and index to access elements/shape functions
	//auto lv = basis.localView();

	//lv.bind(*it);
	
	// Access to shape functions on reference element
	//auto& localBasis = lv.tree().finiteElement().localBasis();
	//auto& localFiniteElement = lv.tree().finiteElement();

	//auto myfunc = [](const auto& x){ return 0.01;};
	
	//std::vector<Dune::FieldVector<double,1>> interpolates;
	//localFiniteElement.localInterpolation(myfunc,interpolates);
	//std::cout << interpolates << '\n';


	std::ofstream file;
	file.open("hermite");
	HermiteFiniteElement<double> localBasis;
	std::vector<Dune::FieldVector<double,1>> Nx;
	for(double x=0; x<=1; x+=0.01)
	{
		localBasis.evaluateFunction(x,Nx);
		file << x;
		for(int i=0;i<localBasis.size();i++)
		{
			file << '\t' << Nx[i];
		}
		file << '\n';
	}
	file.close();
	
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}


