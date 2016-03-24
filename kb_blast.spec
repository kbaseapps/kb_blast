/*
** A KBase module: kb_blast
**
** This module contains 6 methods from BLAST+: BLASTn, BLASTp, BLASTx, tBLASTx, tBLASTn, and PSI-BLAST
** 
*/

module kb_blast {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /* BLAST Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	sequence       input_one_sequence;
	data_obj_name  input_one_name;
	data_obj_name  input_many_name;
        data_obj_name  output_filtered_name;

	float ident_thresh;
	float e_value;
	float bitscore;
	float overlap_fraction;
	float maxaccepts;
	float rounds;  /* for PSI-BLAST */
    } BLAST_Params;


    /* BLAST Output
    */
    typedef structure {
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
    } BLAST_Output;
	

    /*  Methods for BLAST of various flavors of one sequence against many sequences 
    **
    **    overloading as follows:
    **        input_one_id: SingleEndLibrary, Feature, FeatureSet
    **        input_many_id: SingleEndLibrary, FeatureSet, Genome, GenomeSet
    **        output_id: SingleEndLibrary (if input_many is SELib), (else) FeatureSet
    */
    funcdef BLASTn_Search (BLAST_Params params)  returns (BLAST_Output) authentication required;
    funcdef BLASTp_Search (BLAST_Params params)  returns (BLAST_Output) authentication required;
    funcdef BLASTx_Search (BLAST_Params params)  returns (BLAST_Output) authentication required;
    funcdef tBLASTn_Search (BLAST_Params params)  returns (BLAST_Output) authentication required;
    funcdef tBLASTx_Search (BLAST_Params params)  returns (BLAST_Output) authentication required;
};
