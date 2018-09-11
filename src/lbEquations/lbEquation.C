#include <lbEquation.H>


/** Default constructor */

lbEquation::lbEquation( const std::string& name,
			const latticeMesh& mesh_,
			timeOptions& Time_,
			pdfField& pdf_) : ename (name),
					  mesh(mesh_),
					  Time(Time_),
					  _q( mesh.lmodel()->q() ),					  
					  pdf(pdf_){}


/** Default destructor */

lbEquation::~lbEquation() {}



/** Collision process */

const void lbEquation::collision() {}


/** Streamming process */

const void lbEquation::streamming() {}
