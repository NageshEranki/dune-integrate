/*

Extra functions/classes

*/

template<class T>
void printqr(T& quad)
{
	std::cout << "Printing quadrature rules:" << '\n';
	for(int i=0;i<quad.size();i++)
	{
		std::cout << "x" << i << " = " << quad[i].position() << " , w" << i << " = " << quad[i].weight() << '\n';
	}
}

template<class V, class E>
void printsfeval(V& localView, E& element)
{
	localView.bind(*element);
	auto& localBasis = localView.tree().finiteElement().localBasis();
	std::cout << "Evaluating shape functions on Element nodes:\n";
	std::vector<Dune::FieldVector<double,1>> sfevals;
	auto geom = element->geometry();
	for(int i=0;i<geom.corners();i++)
	{
		localBasis.evaluateFunction(geom.local(geom.corner(i)),sfevals);
		for(int j=0;j<localBasis.size();j++)
		{
			std::cout << "d(" << j << "," << i << ") = " << sfevals[j] << '\n';
		}
	}

}


void useless()
{
	std::cout << "THis function does nothing\n";
}	
