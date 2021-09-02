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
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/common/fvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

template<class CoordinateType>
bool isDirichlet(const CoordinateType& xy)
{
	if((xy[0] == 0))// || (xy[0] == 1)) // || (xy[1] == 1))
		return true;
	return false;
}

template<class T>
double dirichletFunc(const T& x)
{
	if(x[0] == 0)
		return 1;
	/*	
	if(x[0] == 1)
		return 0;
	*/
	return 0;
}

template<class T>
double neumannFunc(const T& x)
{
	if(x[0] == 1) 
		return -1;
	return 0;
}



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
	GridType grid({1.0},{3});
	//GridType grid({1.0,1.0},{3,3});

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

	lv.bind(*it);

	auto& localBasis = lv.tree().finiteElement().localBasis();

	std::vector<Dune::FieldMatrix<double,1,dim>> refGrad;

	MatrixType stiffnessMatrix(basis.dimension(),basis.dimension(),50,0.2,Dune::BCRSMatrix<double>::implicit);

	//Begin assembly of Stiffness Matrix
	
	for(const auto& element: elements(gv))
	{
		// Bind Local View to current elem
		lv.bind(element);
		
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
					stiffnessMatrix.entry(lv.index(i),lv.index(j)) 
						+= weight*(grad[i]*grad[j])*geom.integrationElement(pos);
					
				}
			}
		}
	}

	stiffnessMatrix.compress();
	
	//End assembly


	// Load vector init
	VectorType b(basis.dimension());
	b = 0;

	const auto& indexSet = gv.indexSet();
	

	// Assemble contributions from Neumann condition to load vector

	std::vector<Dune::FieldVector<double,1>>  evals;

	const auto& surfquad = Dune::QuadratureRules<double,dim-1>::rule(gv.ibegin(*gv.begin<0>())->geometry().type(),2);

	for(const auto& element: elements(gv))
	{
		if(element.hasBoundaryIntersections())
		{
			lv.bind(element);
			for(auto it=gv.ibegin(element);it!=gv.iend(element);it++)
			{
				auto geom = it->geometry();
				for(const auto& quad: surfquad)
				{
					auto quadPt = quad.position();
					auto weight = quad.weight();
					localBasis.evaluateFunction(element.geometry().local(geom.global(quadPt)),evals);
					for(int i=0; i<localBasis.size(); i++)
					{
						b[lv.index(i)] += weight*neumannFunc(element.geometry().local(geom.global(quadPt)))*geom.integrationElement(quadPt);
					}
					//res += weight*n*dirichletFunc(geom.global(quadPt))*geom.integrationElement(quadPt);
				}
			}
		}

	}
	//End Neumann


	// Visit each vertex to search for dirichlet nodes
	for(const auto& vertex: vertices(gv))
	{
		auto pos = vertex.geometry().corner(0);
		if(isDirichlet(pos))
		{
			// If current node is dirichlet, note down the dof in the solution vector x
			auto dof = indexSet.index(vertex);
			auto* rowptr = &stiffnessMatrix[dof];

			// Replace corresponding row in stiffness matrix with
			// row from identity matrix
			for(auto ptr = rowptr->begin(); ptr != rowptr->end(); ptr++)
				*ptr = (dof == ptr.index()) ? 1 : 0;

			// Add entry to load vector
			b[dof] = dirichletFunc(pos);
		}
	}
	//End Dirichlet


	Dune::printmatrix(std::cout,stiffnessMatrix,"Stiffness Matrix","row");
	Dune::printvector(std::cout,b,"Load Vector","Row");

	Dune::MatrixAdapter<MatrixType,VectorType,VectorType> lop(stiffnessMatrix);
	Dune::SeqILU<MatrixType,VectorType,VectorType> pre(stiffnessMatrix,1.0);
	Dune::CGSolver<VectorType> cg(lop,pre,1e-8,1000,3);
	Dune::InverseOperatorResult stat;	
	VectorType x;
	x = b;
	x=0;
	Dune::printvector(std::cout,x,"Pre-Solution","");
	
	cg.apply(x,b,stat);

	Dune::printvector(std::cout,x,"Solution","");
	
	Dune::VTKWriter<GridView> writer(gv);
	auto func = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis,x);
	writer.addVertexData(func,Dune::VTK::FieldInfo("Temp",Dune::VTK::FieldInfo::Type::scalar,1));

	writer.write("sol-laplace");

	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
