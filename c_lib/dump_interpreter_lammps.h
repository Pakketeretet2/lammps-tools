#ifndef DUMP_INTERPRETER_LAMMPS_H
#define DUMP_INTERPRETER_LAMMPS_H

#include "block_data.h"
#include "dump_reader.h"
#include "text_readers.h"


class dump_interpreter_lammps : public dump_interpreter
{
public:
	dump_interpreter_lammps( const std::string &dname, int file_type = 0 );
	virtual ~dump_interpreter_lammps();
	
	virtual int next_block( block_data &block );

	virtual int next_block_meta( block_data &block );
	virtual int next_block_body( block_data &block );

	virtual bool eof() const
	{
		if ( r ) return r->eof();
		else     return false;
	}
	virtual bool good() const
	{
		if( r ) return r->good();
		else    return false;
	}

	enum dump_styles {
		ATOMIC = 0,
		CUSTOM
	};
	
private:
	enum { BIT_X = 1,
	       BIT_Y = 2,
	       BIT_Z = 3 };

	text_reader *r;
	
	int id_idx;
	int type_idx;
	int x_idx;
	int y_idx;
	int z_idx;
	int mol_idx;

	int n_cols;
	std::vector<std::string> headers;
	
	std::vector<int>         other_cols;
	std::vector<std::string> other_col_headers;

	int atom_style;
	int scaled;
	std::string last_line;
	block_data last_meta;

	int dump_style;

	void set_headers( const std::string &h_line );
};






#endif // DUMP_INTERPRETER_LAMMPS_H
