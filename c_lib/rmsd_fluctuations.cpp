#include "rmsd_fluctuations.h"
#include "id_map.h"
#include "util.h"

#include <algorithm>

using vec_arr_4 = std::vector<std::array<double,4> >;


py_float rmsd_fluctuations_impl( dump_reader &r, std::vector<std::array<double,4> > &rmsd_flucts,
                                 const std::list<long int> *ids, std::list<double> &rmsd_time )
{
	double rmsd_total = 0;
	block_data block;
	rmsd_time.clear();
	
	r.next_block( block );

	rmsd_flucts.resize( block.N );


	// Sort by ids:
	std::vector<std::array<double, 3> > x_avg( block.N );
	id_map im0( block.ids ); // This gives you the indices at time step 0.
	
	std::vector<std::array<double, 3> > x_sor( block.N );
	double alpha = 0.95;

	x_avg = block.x;

	do {
		std::cerr << "Block has " << block.N << " atoms, t = "
		          << block.tstep << "...\n";
		
		if( ids ){
			block = filter_block( block, *ids );
		}
		
		std::cerr << "Done filtering, now got " << block.N
		          << " atoms...\n";
		// Sort the positions into the same order as im0.
		for( int i = 0; i < block.N; ++i ){
			int id  = block.ids[i];
			int idx = im0.id_to_index( id );

			x_sor[idx] = block.x[i];
		}

		// Now perform the moving average:
		double rmsd_curr = 0.0;
		for( int i = 0; i < block.N; ++i ){

			// Determine per-atom rmsd from average:
			double dx = x_avg[i][0] - x_sor[i][0];
			double dy = x_avg[i][1] - x_sor[i][1];
			double dz = x_avg[i][2] - x_sor[i][2];
			double x2 = dx*dx;
			double y2 = dy*dy;
			double z2 = dz*dz;
			
			rmsd_flucts[i][0] = x2;
			rmsd_flucts[i][1] = y2;
			rmsd_flucts[i][2] = z2;
			rmsd_flucts[i][3] = x2 + y2 + z2;

			rmsd_curr += rmsd_flucts[i][3];

			x_avg[i][0] = (1-alpha)*x_avg[i][0] + alpha*x_sor[i][0];
			x_avg[i][1] = (1-alpha)*x_avg[i][1] + alpha*x_sor[i][1];
			x_avg[i][2] = (1-alpha)*x_avg[i][2] + alpha*x_sor[i][2];


			
		}
		std::cerr << "At t = " << block.tstep << ", rmsd = "
		          << rmsd_curr << "\n";

		rmsd_time.push_back( rmsd_curr );
		rmsd_total += rmsd_curr;

	}while( r.next_block( block ) );
	

	return rmsd_total;
}


py_float get_msd( const block_data &b1, const block_data &b0, py_float *msd )
{
	double tmsd = 0.0;
	id_map im1( b1.ids );
	
	for( int i = 0; i < b0.N; ++i ){
		long int id0 = b0.ids[i];
		const std::array<double,3> &x0 = b0.x[i];
		const std::array<double,3> &x1 = b1.x[ im1[id0] ];

		double dx = x0[0] - x1[0];
		double dy = x0[1] - x1[1];
		double dz = x0[2] - x1[2];
		double r2 = dx*dx + dy*dy + dz*dz;

		tmsd += r2;
		msd[i] = r2;
	}
	
	return tmsd;
}
