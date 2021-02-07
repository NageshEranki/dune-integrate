// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

//Grid includes
#include <dune/grid/yaspgrid.hh>
//#include <dune/grid/uggrid.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dune/common/fvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/bcrsmatrix.hh>

template<class T>
double func(T x)
{
	return 1;
}

int main(int argc, char** argv)
{
  try
  {
	//Start!
	
	// Dimension of domain.
	const int dim = 2;

	typedef Dune::FieldMatrix<double,1,1> FType;
	
	// Structured grid from YaspGrid
	typedef Dune::YaspGrid<dim> GridType;
	GridType grid({1.0,1.0},{5,5});

	//Obtain view into grid
	typedef GridType::LeafGridView GridView;
	GridView gv = grid.leafGridView();

	//Pointer to first element in the grid
	auto it = gv.begin<0>();

	//Obtain Quadrature rules from geometry type
	const auto& quad = Dune::QuadratureRules<double,dim>::rule(it->geometry().type(),2);

	typedef Dune::Functions::LagrangeBasis<GridView,1> BasisType;
	BasisType basis(gv);
	
	// Create local views and index to access elements/shape functions
	auto lv = basis.localView();
	auto localIndexSet = basis.localIndexSet();

	lv.bind(*it);
	localIndexSet.bind(lv);
	auto& localBasis = lv.tree().finiteElement().localBasis();

	std::vector<Dune::FieldMatrix<double,1,dim>> refGrad;
	Dune::FieldMatrix<double,4,4> elementMatrix(0);

	Dune::BCRSMatrix<FType> stiffnessMatrix(basis.dimension(),basis.dimension(),6,0.2,Dune::BCRSMatrix<FType>::implicit);

	for(const auto& element: elements(gv))
	{
		// Bind Local View to current elem
		lv.bind(element);
		
		// Bind LocalIndexSet to updated local view.
		// Maps shape functions on ref elem to global dofs
		localIndexSet.bind(lv);
		
		// Geometry instance to evaluate jacobian, determinants, etc
		auto geom = element.geometry();

		
		for(int i=0; i<localBasis.size(); i++)
		{
			for(int j=0; j<localBasis.size(); j++)
			{
				for(const auto& quadPt: quad)
				{
					const auto& pos = quadPt.position();
					const auto& weight = quadPt.weight();

					// Evaluate gradients on ref elem
					localBasis.evaluateJacobian(pos,refGrad);

					std::vector<Dune::FieldVector<double,dim>> grad(refGrad.size());

					auto jacobian = geom.jacobianInverseTransposed(pos);
					
					// Transform gradients on ref elem to current elem
					for(int k=0;k<refGrad.size();k++)
						jacobian.mv(refGrad[k][0],grad[k]);

					// Add Entry to sparse matrix
					stiffnessMatrix.entry(localIndexSet.index(i),localIndexSet.index(j)) 
						+= weight*(grad[i]*grad[j])*geom.integrationElement(pos);
					
				}
			}
		}
	}

	stiffnessMatrix.compress();
	
	Dune::printmatrix(std::cout, stiffnessMatrix, "K","");

	const auto& indexSet = gv.indexSet();
	int count = 0;
	for(const auto& vertex: vertices(gv))
	{
		auto pos = vertex.geometry().corner(0);
		if((pos[1] == 0) || (pos[0] == 0) || (pos[0]==1) || (pos[1]==1))
		{
			std::cout << count 
				  << "------>"
				  << "Dirichlet dof " 
				  << indexSet.index(vertex) 
				  << " @ " <<  vertex.geometry().corner(0) 
				  << std::endl;
			count++;
		}
	}
	std::cout << "Total no of dofs " << gv.size(2) << std::endl;
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
