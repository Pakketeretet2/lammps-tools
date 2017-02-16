#include "neighborize.h"

#ifndef HAVE_LIB_CGAL
// Implementations are no-ops, just to make it compile.
void neighborize_delaunay( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           list *neighs, py_int itype, py_int jtype ) {}

void neighborize_conv_hull( const arr3f &x, py_int N, const arr1i &ids,
                            const arr1i &types, py_int periodic,
                            const py_float *xlo, const py_float *xhi, py_int dims,
                            list *neighs, py_int itype, py_int jtype ) {}

#else
// Actual implementation using CGAL. Requires many includes...
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <fstream>

#include "domain.h"
#include "id_map.h"
#include "my_output.hpp"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

static my_ostream my_out( std::cout );

inline bool in_list( py_int jd, const list &ni )
{
	return std::find( ni.begin(), ni.end(), jd ) != ni.end();
}



void neighborize_delaunay_2d( const arr3f &x, py_int N, const arr1i &ids,
                              const arr1i &types, const py_float *xlo,
                              const py_float *xhi,
                              list *neighs, py_int itype, py_int jtype )
{
	typedef CGAL::Triangulation_vertex_base_with_info_2<py_int, K> Vb;
	typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
	typedef CGAL::Delaunay_triangulation_2<K, Tds> Dt;
	typedef Dt::Point point;

	my_out << "Triangulating the following data set:\n";
	my_out << "i id type x y z\n";
	for( py_int i = 0; i < N; ++i ){
		my_out << i << " " << ids[i] << " " << types[i] << " "
		          << x[i][0] << " " << x[i][1] << " " << x[i][2] << "\n";
	}

	std::vector<std::pair<point, py_int> > pts;
	pts.reserve(N);
	for( py_int i = 0; i < N; ++i ){
		point pt( x[i][0], x[i][1] );
		pts.push_back( std::make_pair( pt, ids[i]) );
	}
	id_map im(ids);

	Dt t;
	t.insert( pts.begin(), pts.end() );

	Dt::Finite_faces_iterator fi = t.faces_begin();
	for( ; fi != t.faces_end(); ++fi ){
		// For each face, find the points in it and add these
		// to each others neigh lists...
		py_int id = fi->vertex(0)->info();
		py_int jd = fi->vertex(1)->info();
		py_int kd = fi->vertex(2)->info();

		my_out << "Apparently ids " << id << ", " << jd << " and "
		          << kd << " share a face. Their positions are:\n";
		py_int i = im[id], j = im[jd], k = im[kd];
		
		my_out << "    ( " << x[i][0] << ", " << x[i][1] << " ),\n"
		          << "    ( " << x[j][0] << ", " << x[j][1] << " ),\n"
		          << "    ( " << x[k][0] << ", " << x[k][1] << " ),\n";
		if( !in_list( j, neighs[i] ) ) neighs[i].push_back(j);
		if( !in_list( k, neighs[i] ) ) neighs[i].push_back(k);

		if( !in_list( i, neighs[j] ) ) neighs[j].push_back(j);
		if( !in_list( k, neighs[j] ) ) neighs[j].push_back(k);
		
		if( !in_list( i, neighs[k] ) ) neighs[k].push_back(i);
		if( !in_list( j, neighs[k] ) ) neighs[k].push_back(j);
	}
	// Print a gnuplottable test file:
	std::ofstream out( "test_2d_del_gnuplot_lines.dat" );
	std::ofstream out2( "test_2d_del_gnuplot_points.dat" );
	for( py_int i = 0; i < N; ++i ){
		for( py_int j : neighs[i] ){
			out << x[i][0] << " " << x[i][1] << "\n"
			    << x[j][0] << " " << x[j][1] << "\n\n";
		}
		out2 << x[i][0] << " " << x[i][1] << "\n";
	}
}

void neighborize_delaunay_3d( const arr3f &x, py_int N, const arr1i &ids,
                              const arr1i &types, const py_float *xlo,
                              const py_float *xhi,
                              list *neighs, py_int itype, py_int jtype )
{}



void neighborize_delaunay_p_2d( const arr3f &x, py_int N, const arr1i &ids,
                                const arr1i &types, const py_float *xlo,
                                const py_float *xhi,
                                py_int periodic, list *neighs, py_int itype, py_int jtype )
{
	typedef CGAL::Periodic_2_triangulation_filtered_traits_2<K> Gt;
	typedef CGAL::Periodic_2_triangulation_vertex_base_2<Gt> Vb_base;
	typedef CGAL::Triangulation_vertex_base_with_info_2<py_int, Gt, Vb_base> Vb;
	typedef CGAL::Periodic_2_triangulation_face_base_2<Gt> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt, Tds> PDt;
	typedef PDt::Point point;

	typedef PDt::Iso_rectangle rect;
	rect domain( xlo[0], xlo[1], xhi[0], xhi[1] );

	
	std::vector<std::pair<point, py_int> > pts(N);
	for( py_int i = 0; i < N; ++i ){
		const py_float *xi = x[i];
		pts.push_back( std::make_pair( point(xi[0], xi[1]), ids[i] ) );
		my_out << "Point " << i << " = ( " << xi[0] << ", "
		          << xi[1] << " ).\n";
		if( xi[0] < xlo[0] || xi[0] > xhi[0] ||
		    xi[1] < xlo[1] || xi[1] > xhi[1] ){
			std::cerr << "Point xi outside of domain! Domain = [ "
			          << xlo[0] << ", " << xhi[0] << " ] x [ "
			          << xlo[1] << ", " << xhi[1] << " ]\n";
		}

	}

	PDt t( pts.begin(), pts.end(), domain );

	id_map im( ids );

	PDt::Finite_faces_iterator fi = t.faces_begin();
	for( ; fi != t.faces_end(); ++fi ){	
		// For each face, find the points in it and add these
		// to each others neigh lists...
		py_int id = fi->vertex(0)->info();
		py_int jd = fi->vertex(1)->info();
		py_int kd = fi->vertex(2)->info();

		
		// check if points are ok (you can abuse the fact that all ids are >0,
		// while the py_ints are default-constructed to 0).
		bool id_no, jd_no, kd_no;
		id_no = !id;
		jd_no = !jd;
		kd_no = !kd;
		int i = im[id], j = im[jd], k = im[kd];
		
		if( !id_no && !jd_no ){
			if( !in_list( j, neighs[i] ) ) neighs[i].push_back(j);	
			if( !in_list( i, neighs[j] ) ) neighs[j].push_back(i);
		}

		if( !jd_no && !kd_no ){
			if( !in_list( k, neighs[j] ) ) neighs[j].push_back(k);
			if( !in_list( j, neighs[k] ) ) neighs[k].push_back(j);
		}

		if( !id_no && !kd_no ){
			if( !in_list( k, neighs[i] ) ) neighs[i].push_back(k);
			if( !in_list( i, neighs[k] ) ) neighs[k].push_back(i);
		}


	}
	// Print a gnuplottable test file:
	std::ofstream out( "test_2d_p_del_gnuplot_lines.dat" );
	std::ofstream out2( "test_2d_p_del_gnuplot_points.dat" );
	for( py_int i = 0; i < N; ++i ){
		for( py_int j : neighs[i] ){
			out << x[i][0] << " " << x[i][1] << "\n"
			    << x[j][0] << " " << x[j][1] << "\n\n\n";
		}
		out2 << x[i][0] << " " << x[i][1] << "\n";
	}
}


void neighborize_delaunay_p_3d( const arr3f &x, py_int N, const arr1i &ids,
                                const arr1i &types, const py_float *xlo, const py_float *xhi,
                                py_int periodic, list *neighs, py_int itype, py_int jtype )
{

}



void neighborize_delaunay( const arr3f &x, py_int N, const arr1i &ids,
                           const arr1i &types, py_int periodic,
                           const py_float *xlo, const py_float *xhi, py_int dims,
                           list *neighs, py_int itype, py_int jtype )
{
	if( (periodic != PERIODIC_NONE) && (periodic != PERIODIC_FULL) ){
		std::cerr << "Cannot use Delaunay for semi-periodic systems!\n";
		return;
	}
	if( dims == 2 ){
		if( periodic ){
			my_out << "Triangulating (2D periodic Delaunay)...\n";
			neighborize_delaunay_p_2d(x,N,ids,types,xlo,xhi, periodic, neighs,
			                          itype, jtype );
		}else{
			my_out << "Triangulating (2D nonperiodic Delaunay)...\n";
			neighborize_delaunay_2d(x,N,ids,types,xlo,xhi,neighs,
			                        itype, jtype );
		}
	}else{
		if( periodic ){
			my_out << "Triangulating (3D periodic Delaunay)...\n";
			neighborize_delaunay_p_3d(x,N,ids,types,xlo,xhi, periodic, neighs,
			                          itype, jtype );
		}else{
			my_out << "Triangulating (3D nonperiodic Delaunay)...\n";
			neighborize_delaunay_3d(x,N,ids,types,xlo,xhi,neighs,
			                        itype, jtype );
		}
	}
}


void neighborize_conv_hull( const arr3f &x, py_int N, const arr1i &ids,
                            const arr1i &types, py_int periodic,
                            const py_float *xlo, const py_float *xhi, py_int dims,
                            list *neighs, py_int itype, py_int jtype )
{
	typedef CGAL::Polyhedron_3<K> poly;
	typedef K::Point_3 point;

	std::map<point, py_int> pt_to_idx;

	// Some checks on points:
	double R_avg = 0.0;
	double R_var = 0.0;
	double nom = 1.0/N;

	py_int max_id = 0;
	
	for( py_int i = 0; i < N; ++i ){
		double r2 = x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2];
		R_avg += sqrt(r2) * nom;
		max_id = ( max_id < ids[i] ) ? ids[i] : max_id;
	}
	nom = 1.0 / ( N - 1 );
	for( py_int i = 0; i < N; ++i ){
		const py_float *xi = x[i];
		double r = sqrt( xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2] );
		double dr = r - R_avg;
		R_var += dr*dr *nom;
	}
	double R_std = sqrt(R_var);
	const double tol = 1e-4; // Because of poor precision in text files.
	if( R_std > tol ){
		std::cerr << "Warning! Standard deviation in radius = " << R_std
		          << " > tol " << tol << ". Average = " << R_avg
		          << ".\n";
	}

	std::vector<point> pts;
	pts.reserve(N);
	for( py_int i = 0; i < N; ++i ){
		point p( x[i][0], x[i][1], x[i][2] );
		pts.push_back(p);
		pt_to_idx[p] = i;
	}

	std::ofstream out ( "test_conv_hull_gnuplot_lines.dat"  );
	std::ofstream out2( "test_conv_hull_gnuplot_points.dat" );
	
	
	poly hull;
	CGAL::convex_hull_3( pts.begin(), pts.end(), hull );
	poly::Halfedge_iterator ei = hull.halfedges_begin();
	for( ; ei != hull.halfedges_end(); ++ei ){
		// For each facet, find the points in it and add these
		// to each others neigh lists...
		point p1 = ei->vertex()->point();
		if( ei->next() == hull.halfedges_end() ) break;
		
		point p2 = ei->next()->vertex()->point();
		py_int i = pt_to_idx[ p1 ];
		py_int j = pt_to_idx[ p2 ];
		py_int id1 = ids[i];
		py_int id2 = ids[j];

		if( !in_list( j, neighs[i] ) ) neighs[i].push_back(j);
		if( !in_list( i, neighs[j] ) ) neighs[j].push_back(i);

		out << p1[0] << " " << p1[1] << " " << p1[2] << "\n"
		    << p2[0] << " " << p2[1] << " " << p2[2] << "\n\n\n";
		out2 << p1[0] << " " << p1[1] << " " << p1[2] << "\n";
	}
}


#endif // HAVE_LIB_CGAL
