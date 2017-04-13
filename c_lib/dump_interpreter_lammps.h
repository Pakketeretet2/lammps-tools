#ifndef DUMP_INTERPRETER_LAMMPS_H
#define DUMP_INTERPRETER_LAMMPS_H

#include "block_data.h"
#include "dump_reader.h"

class dump_interpreter_lammps : public dump_interpreter
{
public:
	dump_interpreter_lammps() : id_idx(-1), type_idx(-1), x_idx(-1),
	                            y_idx(-1), z_idx(-1), mol_idx(-1),
	                            n_cols(0), atom_style( atom_styles::ATOMIC ),
	                            scaled(0)
	{}
	virtual ~dump_interpreter_lammps(){}
	virtual int next_block( reader_core *r, block_data &block );

	virtual int next_block_meta( reader_core *r, block_data &block );
	virtual int next_block_body( reader_core *r, block_data &block );
	
	
private:
	enum { BIT_X = 1,
	       BIT_Y = 2,
	       BIT_Z = 3 };
	
	int id_idx;
	int type_idx;
	int x_idx;
	int y_idx;
	int z_idx;
	int mol_idx;

	int n_cols;
	std::vector<std::string>          headers;
	
	std::vector<int>                  other_cols;
	std::vector<std::string>          other_col_headers;

	int atom_style;

	int scaled;

	void set_headers( const std::string &h_line );
};






#endif // DUMP_INTERPRETER_LAMMPS_H
