#include <energyNormalHeatFlux.H>

using namespace std;


/** Constructor */

energyNormalHeatFlux::energyNormalHeatFlux( const std::string& eqName,
				    const std::string& bdName,
				    const latticeMesh& mesh,
				    const scalarField& rho,
				    const scalarField& T,
				    const vectorField& U,
				    pdfField& pdf )

    : energyFixedGradT( eqName, bdName, mesh, rho, T, U, pdf ) {


    // Temperature limit

    _Tlimit.resize( _nodes.size() );

    dictionary dict("start/boundaries");    

    scalar tl = dict.lookUpOrDefault<scalar>( eqName + "/" + bdName + "/limit", 1000 );

    for( uint i = 0 ; i < _bndVal.size() ; i++ )
	_Tlimit[i] = tl;    


}




/** Destructor */

energyNormalHeatFlux::~energyNormalHeatFlux() {}



/** Update pdf field */

void energyNormalHeatFlux::update( const energyEquation* eeq ) {


    // First compute value over boundary according to _grad

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	
    	// _bndVal[i] = _T.at(_nodes[i]) + _grad / eeq->thermalCond(_nodes[i]);
    	// _bndVal[i] = _T.at(_nbid[i]) + _grad / eeq->thermalCond(_nbid[i]);
    	_bndVal[i] = ( 4*_T.at(_nbid[i]) - _T.at(_snbid[i]) + 2*_grad[i] / eeq->thermalCond(_nbid[i]) ) / 3;
	// cout << _nodes[i] << " " << _nbid[i] << " " << _snbid[i] << "     " << _T.at(_nbid[i]) << " " << _T.at(_snbid[i]) << " " << _bndVal[i] << endl;


	if( _bndVal[i] > _Tlimit[i] )
	    _bndVal[i] = _Tlimit[i];
	
    }
    

    
    energyFixedT::update( eeq );

}
