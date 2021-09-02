// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
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

template<class T>
void printqr(const T& quad);

template<class V, class E>
void printsfeval(V& localView, const E& element);

template<class V, class E>
void printsfjac(V& localView, const E& element);

template<class E>
void printgeomprop(const E& element);

void useless();

template<class BasisType, class ElementType, class MatrixType,  class Quadrature>
void buildElementMatrix(BasisType& basis, ElementType& element, const Quadrature& quadrule, MatrixType& elemK);

template<class CoordinateType>
double f(const CoordinateType& x);

template<class BasisType, class ElementType, class VectorType,  class Quadrature>
void buildElementLoadVector1(BasisType& basis, ElementType& element, const Quadrature& quadrule, VectorType& elemF);

template<class GridViewType, class MatrixType, class VectorType>
void assembleDirichlet(const GridViewType& gv, MatrixType& K , VectorType& F );

template<class CoordinateType>
bool isDirichlet(const CoordinateType& x);

template<class CoordinateType>
double dirichletFunc(const CoordinateType& x);

template<class T>
double neumannFunc(const T& x);

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
	printqr(quad);

	typedef Dune::Functions::LagrangeBasis<GridView,1> BasisType;
	BasisType basis(gv);
	
	// Create local views and index to access elements/shape functions
	auto lv = basis.localView();

	lv.bind(*it);
	
	// Access to shape functions on reference element
	auto& localBasis = lv.tree().finiteElement().localBasis();

	// Sanity check
	//printsfeval(lv,*it);
	//printsfjac(lv,*it);
	//printgeomprop(*it);

	// Create global stiffness matrix & global force vector
	MatrixType K(basis.size(),basis.size(),3,0.2,MatrixType::implicit);
	VectorType F(basis.size());

	// Create container for element stiffness matrix & element force vector
	const int test = localBasis.size();
	std::cout << "test = " << test << "\n";
	Dune::FieldMatrix<double,2*dim,2*dim> elemK(0);
	Dune::FieldVector<double,2*dim> elemF(0);

	for(const auto element: elements(gv))
	{
		// Iterate over elements and compute
		// respective element stiffness matrices
		// and force vectors

		elemK = {0};
		elemF = {0};
		buildElementLoadVector1(basis,element,quad,elemF);
		buildElementMatrix(basis,element,quad,elemK);
		lv.bind(element);
		for(int i=0;i<localBasis.size();i++)
		{
			for(int j=0;j<localBasis.size();j++)
			{
				// Fit element stiffness matrix
				// into global stiffness matrix
				K.entry(lv.index(i),lv.index(j)) += elemK[i][j];
			}

			// Fit element force vector
			// into global force vector
			F[lv.index(i)] += elemF[i];
		}
	}

	K.compress();

	//Assemble dirichlet condition here
	assembleDirichlet(gv, K , F );

	// Experimental!
	F[0] += 1;
	//
	

	Dune::printmatrix(std::cout,K,"Stiffness Matrix","");
	Dune::printvector(std::cout,F,"Load Vector","");

	// Solver config
	Dune::MatrixAdapter<MatrixType,VectorType,VectorType> lop(K);
	Dune::SeqILU<MatrixType,VectorType,VectorType> pre(K,1.0);
	Dune::CGSolver<VectorType> cg(lop,pre,1e-8,1000,3);
	Dune::InverseOperatorResult stat;	
	VectorType x(basis.size());
	x = F;
	x=0;
	
	// Solve!
	cg.apply(x,F,stat);

	Dune::printvector(std::cout,x,"Solution","");

	// Write one dimensional solution
	// to file, handle with Gnuplot
	Dune::GnuplotWriter writer(gv);
	writer.addVertexData(x,"Nodal solution");
	writer.write("1d_solution");

	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

template<class BasisType, class ElementType, class MatrixType,  class Quadrature>
void buildElementMatrix(BasisType& basis, ElementType& element, const Quadrature& quadrule, MatrixType& elemK )
{
	auto lv = basis.localView();
	lv.bind(element);
	auto& localBasis = lv.tree().finiteElement().localBasis();
	std::vector<Dune::FieldMatrix<double,1,ElementType::dimension>> refGrad;
	for(const auto& quad: quadrule)
	{
		auto geom = element.geometry();
		auto pos = quad.position();
		auto weight = quad.weight();
		localBasis.evaluateJacobian(pos,refGrad);
		std::vector<Dune::FieldVector<double,ElementType::dimension>> grad;
		grad.resize(refGrad.size());
		auto jac = geom.jacobianInverseTransposed(pos);
		for(int k=0;k<refGrad.size();k++)
			jac.mv(refGrad[k],grad[k]);
		for(int i=0;i<localBasis.size();i++)
		{
			for(int j=0;j<localBasis.size();j++)
			{
				elemK[i][j] += weight*grad[i]*grad[j]*geom.integrationElement(pos);
			}
		}
	}
}

template<class BasisType, class ElementType, class VectorType,  class Quadrature>
void buildElementLoadVector1(BasisType& basis, ElementType& element, const Quadrature& quadrule, VectorType& elemF )
{
	/// Contribution of force term
	/// To element load vector
	
	auto lv = basis.localView();
	lv.bind(element);
	auto& localBasis = lv.tree().finiteElement().localBasis();
	auto geom = element.geometry();
	std::vector<Dune::FieldVector<double,1>> sfvals;
	for(const auto& quad: quadrule)
	{
		auto pos = quad.position();
		auto weight = quad.weight();
		localBasis.evaluateFunction(pos,sfvals);
		for(int i=0;i<localBasis.size();i++)
		{
			elemF[i] += weight*sfvals[i]*f(geom.global(pos))*geom.integrationElement(pos);
		}
	}
}

template<class GridViewType, class MatrixType, class VectorType>
void assembleDirichlet(const GridViewType& gv, MatrixType& K , VectorType& F )
{
	auto& indexSet = gv.indexSet();
	for(const auto& vertex: vertices(gv))
	{
		auto pos = vertex.geometry().corner(0);
		if(isDirichlet(pos))
		{
			// If current node is dirichlet, note down the dof in the solution vector x
			auto dof = indexSet.index(vertex);
			auto* rowptr = &K[dof];

			// Replace corresponding row in stiffness matrix with
			// row from identity matrix
			for(auto ptr = rowptr->begin(); ptr != rowptr->end(); ptr++)
				*ptr = (dof == ptr.index()) ? 1 : 0;

			// Add entry to load vector
			F[dof] = dirichletFunc(pos);
		}
	}
}

template<class CoordinateType>
bool isDirichlet(const CoordinateType& x)
{
	if(x[0] == 1)
	{
		return true;
	}
	return false;
}

template<class CoordinateType>
double dirichletFunc(const CoordinateType& x)
{
	return 1;
}

template<class CoordinateType>
double f(const CoordinateType& x)
{
	//return sin(x[0]);
	return x[0];
}

template<class T>
double neumannFunc(const T& x)
{
	if(x[0] == 0) 
		return 1;
	return 0;
}

template<class T>
void printqr(const T& quad)
{
	std::cout << "Printing quadrature rules:" << '\n';
	for(int i=0;i<quad.size();i++)
	{
		std::cout << "x" << i << " = " << quad[i].position() << " , w" << i << " = " << quad[i].weight() << '\n';
	}
	std::cout << "=============================================================\n";
}

template<class V, class E>
void printsfeval(V& localView, const E& element)
{
	localView.bind(element);
	auto& localBasis = localView.tree().finiteElement().localBasis();
	std::cout << "Evaluating shape functions on Element nodes:\n";
	std::vector<Dune::FieldVector<double,E::mydimension>> sfevals;
	auto geom = element.geometry();
	for(int i=0;i<geom.corners();i++)
	{
		localBasis.evaluateFunction(geom.local(geom.corner(i)),sfevals);
		for(int j=0;j<localBasis.size();j++)
		{
			std::cout << "d(" << j << "," << i << ") = " << sfevals[j] << '\n';
		}
	}
	std::cout << "=============================================================\n";

}

template<class V, class E>
void printsfjac(V& localView, const E& element)
{
	localView.bind(element);
	auto& localBasis = localView.tree().finiteElement().localBasis();
	std::vector<Dune::FieldMatrix<double,1,E::mydimension>> sfjac;
	auto geom = element.geometry();
	localBasis.evaluateJacobian(geom.local(geom.center()),sfjac);
	std::cout << "Evaluating shape function derivatives @ element center:" << geom.local(geom.center()) << '\n';
	for(int i=0; i<localBasis.size(); i++)
	{
		Dune::printmatrix(std::cout,sfjac[i],"Shape function derivatives","");
	}
	std::cout << "=============================================================\n";
}

template<class E>
void printgeomprop(const E& element)
{
	auto geom = element.geometry();
	std::cout << "Evaluating Jacobian Inverse Transposed @ " << geom.local(geom.center()) << '\n';
	Dune::printmatrix(std::cout,geom.jacobianInverseTransposed(geom.local(geom.center())),"","");
	std::cout << "=============================================================\n";
	std::cout << "Evaluating integration element @ " << geom.local(geom.center()) << '\n';
	std::cout << geom.integrationElement(geom.local(geom.center())) << '\n';
	std::cout << "=============================================================\n";
}
void useless()
{
	std::cout << "THis function does nothing\n";
}
