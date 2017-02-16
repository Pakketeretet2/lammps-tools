#include "dump_to_povray.h"
#include "dump_reader.h"

#include <fstream>

extern "C" {

void dump_to_povray( const char *dump_fname, const char *povray_fname )
{
	if (!povray_fname){
		povray_fname = "/dev/stdout";
	}
	
	dump_reader d( dump_fname );
	std::ofstream out( povray_fname );
	

	block_data b;
	while( d.next_block(b) ){

		// Calculate a good camera position:
		double x_cam[3];
		x_cam[0] = 4*b.xlo[0];
		x_cam[1] = 4*b.xlo[1];
		x_cam[2] = 4*b.xlo[2];
		
		
		out << "#include \"colors.inc\"\n";
		out << "\n";
		//out << "Width=1600\n"
		//    << "Height=1200\n";
		out << "camera {\n"
		    << "  location <" << x_cam[0] << ", " << x_cam[1] << ", " << x_cam[2] << ">\n"
		    << "  look_at 0\n"
		    << "  angle 36\n"
		    << "}\n"
		    << "light_source { <500,  500, -1000> White }\n"
		    << "light_source { <  0, 1500, -1000> White }\n"
		    << "\n";

		for( int i = 0; i < b.N; ++i ){
			out << "sphere { <0, 0, 0>, 0.5\n"
			    << "  texture { pigment { color White }\n"
				// << "            normal {bumps 0.5 scale 0.05 }\n"
			    << "            finish {ambient 0.1 diffuse 0.9 phong 1 "
			    << "                    irid { 0.25 thickness 0.2 turbulence 0.7} }\n"
			    << "          }\n"
			    << "  translate " << b.x[i][0] << "*x\n"
			    << "  translate " << b.x[i][1] << "*y\n"
			    << "  translate " << b.x[i][2] << "*z\n"
			    << "}\n"
			    << "\n";
		}
		break;
	}
}

}
