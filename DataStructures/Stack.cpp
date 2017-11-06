
#include "Stack.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <math.h>
using namespace NeuroProof;
using namespace std;

#define SIMULATE_SECOND_CHANNEL true

#define BLOCK_SIZE_LIMIT 200
#define BLOCK_SIZE_LIMIT_Z 10

float array_items[256];

void _set_basic_features(FeatureMgr* fman){

    double hist_percentiles[]={0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> percentiles(hist_percentiles, hist_percentiles+sizeof(hist_percentiles)/sizeof(double));		

    // ** for Toufiq's version ** feature_mgr->add_inclusiveness_feature(true);  	
    fman->add_moment_feature(4,true);	
    fman->add_hist_feature(25,percentiles,false); 	

    //cout << "Number of features, channels:" << feature_mgr->get_num_features()<< ","<<feature_mgr->get_num_channels()<<"\n";	

}

Label find_max(Label* data, const size_t* dims){
    Label max=0;  	
    size_t plane_size = dims[1]*dims[2];
    size_t width = dims[2];	 	
    for(size_t i=0;i<dims[0];i++){
	for (size_t j=0;j<dims[1];j++){
	    for (size_t k=0;k<dims[2];k++){
		size_t curr_spot = i*plane_size + j * width + k;
		if (max<data[curr_spot])
		    max = data[curr_spot];	
	    }
	}
    }	
    return max;	
}

void Stack::compute_contingency_table(){

    if(!gtruth)
	return;		

    size_t i,j,k;
    const size_t dimn[]={depth-2*padding,height-2*padding,width-2*padding};
    Label *watershed_data = get_label_volume(); 		
	
    Label ws_max = find_max(watershed_data,dimn);  	

    contingency.clear();	
    //contingency.resize(ws_max+1);
    unsigned long ws_size = dimn[0]*dimn[1]*dimn[2];	
    for(unsigned long itr=0; itr<ws_size ; itr++){
        Label wlabel = watershed_data[itr];
	Label glabel = gtruth[itr];
	multimap<Label, vector<LabelCount> >::iterator mit;
 	mit = contingency.find(wlabel);
	if (mit != contingency.end()){
	    vector<LabelCount>& gt_vec = mit->second;
	    for (j=0; j< gt_vec.size();j++)
	        if (gt_vec[j].lbl == glabel){
		    (gt_vec[j].count)++;
		    break;
	  	}
	    if (j==gt_vec.size()){
	        LabelCount lc(glabel,1);
	        gt_vec.push_back(lc);	
	    }
	}
	else{
	    vector<LabelCount> gt_vec;	
	    gt_vec.push_back(LabelCount(glabel,1));	
	    contingency.insert(make_pair(wlabel, gt_vec));	
	}
    }		

    delete[] watershed_data;	
}
void Stack::compute_vi(){

    if(!gtruth)
	return;		

    int j, k;

    compute_contingency_table();

    double sum_all=0;

    int nn = contingency.size();
	
    multimap<Label, double>  wp;	 		
    multimap<Label, double>  gp;	 		

    for(multimap<Label, vector<LabelCount> >::iterator mit = contingency.begin(); mit != contingency.end(); ++mit){
	Label i = mit->first;
	vector<LabelCount>& gt_vec = mit->second; 	
	wp.insert(make_pair(i,0.0));
	for (j=0; j< gt_vec.size();j++){
	    unsigned int count = gt_vec[j].count;
	    Label gtlabel = gt_vec[j].lbl;
	
	    (wp.find(i)->second) += count;	

	    if(gp.find(gtlabel) != gp.end())
		(gp.find(gtlabel)->second) += count;
	    else
		gp.insert(make_pair(gtlabel,count));	

	    sum_all += count;
	}
	int tt=1;
    }

    double HgivenW=0;
    double HgivenG=0;
     				
    for(multimap<Label, vector<LabelCount> >::iterator mit = contingency.begin(); mit != contingency.end(); ++mit){
	Label i = mit->first;
	vector<LabelCount>& gt_vec = mit->second; 	
	for (j=0; j< gt_vec.size();j++){
	    unsigned int count = gt_vec[j].count;
	    Label gtlabel = gt_vec[j].lbl;
	    
 	    double p_wg = count/sum_all;
	    double p_w = wp.find(i)->second/sum_all; 	
	    
            HgivenW += p_wg* log(p_w/p_wg);

	    double p_g = gp.find(gtlabel)->second/sum_all;	

	    HgivenG += p_wg * log(p_g/p_wg);	

	}
    }

    printf("MergeSplit: (%f, %f)\n",HgivenW, HgivenG); 		
}


void Stack::compute_groundtruth_assignment(){

    if(!gtruth)
	return;		

    size_t j,k;
    const size_t dimn[]={depth-2*padding,height-2*padding,width-2*padding};

    Label gt_max = find_max(gtruth,dimn);  	

    compute_contingency_table();
		
    assignment.clear();  	
    for(multimap<Label, vector<LabelCount> >::iterator mit = contingency.begin(); mit != contingency.end(); ++mit){
	Label i = mit->first;
	vector<LabelCount>& gt_vec = mit->second; 	
	size_t max_count=0;
	Label max_label=0;
	for (j=0; j< gt_vec.size();j++){
	    if (gt_vec[j].count > max_count){
		max_count = gt_vec[j].count;
		max_label = gt_vec[j].lbl;
	    }	
	}
	assignment.insert(make_pair(i,max_label));
    }	


    printf("gt label determined for %d nodes\n", assignment.size());
}
void Stack::modify_assignment_after_merge(Label node_keep, Label node_remove){

    if(!gtruth)
	return;		
  	
    Label node_remove_gtlbl = assignment.find(node_remove)->second;
    size_t node_remove_gtoverlap = 0;
    size_t j;

    multimap<Label, vector<LabelCount> >::iterator mit; 
    mit = contingency.find(node_remove);	
    vector<LabelCount>& gt_vecr = mit->second; 	
    for (j=0; j < gt_vecr.size();j++){
	if (gt_vecr[j].lbl == node_remove_gtlbl){
	    node_remove_gtoverlap = gt_vecr[j].count;
	    break;
	}	
    }
    if (j == gt_vecr.size())
	printf("Something wrong with gt assignment of node remove %d\n",node_remove); 	


    mit = contingency.find(node_keep);	
    vector<LabelCount>& gt_veck = mit->second; 	

    for (j=0; j< gt_veck.size();j++){
	if (gt_veck[j].lbl == node_remove_gtlbl){
	    (gt_veck[j].count) += node_remove_gtoverlap;
	    break;
	}	
    }
    if (j == gt_veck.size()){ // i.e.,  a false merge
	//printf("Node keep %d does not overlap with node remove %d gt \n",node_keep, node_remove); 	
	
	LabelCount lc(node_remove_gtlbl,node_remove_gtoverlap);
	gt_veck.push_back(lc);
    }

    size_t max_count=0;
    Label max_label=0;
    for (j=0; j< gt_veck.size();j++){
	if (gt_veck[j].count > max_count){
	    max_count = gt_veck[j].count;
	    max_label = gt_veck[j].lbl;
	}	
    }
    assignment.find(node_keep)->second = max_label;


}





/*********************************************************************************/




#define TFK_ARRAY_ACCESS(array, stride, x,y,z) array[stride[0]*(x)+stride[1]*(y)+stride[2]*(z)]

void Stack::build_rag_loop(Rag<Label>* rag, FeatureMgr* feature_man, 
			   int x_start, int x_end, int y_start, int y_end, 
			   int z_start, int z_end, bool use_mito_prob) {
    if (feature_man->simulate_channel_1()) {
	Stack::build_rag_loop_tfk(rag, feature_man, x_start, x_end, y_start, y_end,
	    z_start, z_end, use_mito_prob);
    } else {
	Stack::build_rag_loop_leek(rag, feature_man, x_start, x_end, y_start, y_end,
	    z_start, z_end, use_mito_prob);
    }
}

void Stack::build_rag_loop_tfk(Rag<Label>* rag, FeatureMgr* feature_man, 
                                            int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob)
{
    unsigned int maxx = width-1; 
    unsigned int maxy = height-1; 
    unsigned int maxz = depth-1; 

    unsigned int size_estimate = (z_end-z_start)*(y_end-y_start)*(x_end-x_start);
    tfk_edge* sorted_edges = (tfk_edge*) malloc(sizeof(tfk_edge)*size_estimate*7);
    int sorted_edges_len = 0;

    unsigned int * label_vol_array = watershed;
    unsigned int plane_size = width * height;
    uint64_t label_vol_stride[3] = {1,width,plane_size}; 
    tfk_edge edge;

    if (x_start < 0) x_start = 0;
    if (y_start < 0) y_start = 0;
    if (z_start < 0) z_start = 0;

    bool avoid_conditionals = true;
    if (z_start - 1 < 0 || y_start -1 < 0 || x_start - 1 < 0 || z_end+1>maxz || y_end+1 > maxy || x_end+1 > maxx) {
      avoid_conditionals = false;
    }


    unsigned int spots[7];

    if (avoid_conditionals) {
      unsigned int* spots_ptr[7];
      int64_t curr_spot = x_start + y_start*width + z_start*plane_size;
      spots_ptr[0] = watershed + curr_spot;
      spots_ptr[1] = watershed + curr_spot-1;
      spots_ptr[2] = watershed + curr_spot+1;
      spots_ptr[3] = watershed + curr_spot-width;
      spots_ptr[4] = watershed + curr_spot+width;
      spots_ptr[5] = watershed + curr_spot-plane_size;
      spots_ptr[6] = watershed + curr_spot+plane_size;

      for (unsigned int z = z_start; z < z_end; z++) {
          int z_spot = z * plane_size;
          spots_ptr[5] += plane_size;
          spots_ptr[6] += plane_size;
          for (unsigned int y = y_start; y < y_end; y++) {
              int y_spot = y * width;
              spots_ptr[3] += width;
              spots_ptr[4] += width;
              for (unsigned int x = x_start; x < x_end; x++) {
                  curr_spot = x + y_spot + z_spot;
                  spots_ptr[0] += 1;
                  spots_ptr[1] += 1;
                  spots_ptr[2] += 1;
                  for (int i = 0; i < 7; i++) {
                    spots[i] = *(spots_ptr[i]);
                  }


                  unsigned char pred = prediction_array[0][curr_spot];

                  int my_label = spots[0];

                  // sort.
                  for (int j = 0; j < 7; j++) {
                     int insert = spots[j];
                     int hole = j;
                     while (hole > 0 && insert > spots[hole-1]) {
                       spots[hole] = spots[hole-1];
                       hole--;
                     }
                     spots[hole] = insert;
                  }
          
                  int last_label = my_label;
                  bool boundary = false;

                  for (int i = 0; i < 7; i++) {
                    if (spots[i] == last_label || spots[i] == my_label) continue;
                    if (!spots[i]) {
                      boundary = true;
                      last_label = spots[i];
                      continue;
                    }
                    edge.id1 = my_label;
                    edge.id2 = spots[i];
                    edge.pred = pred;
                    edge.boundary = false;
                    sorted_edges[sorted_edges_len++] = edge;
                    last_label = spots[i];
                  }
                  // add a self edge.
                  edge.id1 = my_label;
                  edge.id2 = my_label;
                  edge.pred = pred;
                  edge.boundary = (boundary || !my_label);
                  sorted_edges[sorted_edges_len++] = edge;
              }
          }
      }
    } else {
      for (unsigned int z = z_start; z < z_end; z++) {
          int z_spot = z * plane_size;
          //if (z < 0 || z > maxz) continue;
          for (unsigned int y = y_start; y < y_end; y++) {
              //if (y < 0 || y > maxy) continue;
              int y_spot = y * width;
              for (unsigned int x = x_start; x < x_end; x++) {
                  //if (x < 0 || x > maxx) continue;
                  int64_t curr_spot = x + y_spot + z_spot;
                  spots[0] = watershed[curr_spot];
                  spots[1] = x >= 0 ? watershed[curr_spot-1] : 0;
                  spots[2] = x+1 <= maxx ? watershed[curr_spot+1] : 0;
                  spots[3] = y-1 >=0 ? watershed[curr_spot-width] : 0;
                  spots[4] = y+1 <= maxy ? watershed[curr_spot+width] : 0;
                  spots[5] = z-1 >= 0 ? watershed[curr_spot-plane_size] : 0;
                  spots[6] = z+1 <= maxz ? watershed[curr_spot+plane_size] : 0;

                  unsigned char pred = prediction_array[0][curr_spot];

                  int my_label = spots[0];

                  // sort.
                  for (int j = 0; j < 7; j++) {
                     int insert = spots[j];
                     int hole = j;
                     while (hole > 0 && insert > spots[hole-1]) {
                       spots[hole] = spots[hole-1];
                       hole--;
                     }
                     spots[hole] = insert;
                  }
          
                  int last_label = my_label;
                  bool boundary = false;

                  for (int i = 0; i < 7; i++) {
                    if (spots[i] == last_label || spots[i] == my_label) continue;
                    if (!spots[i]) {
                      boundary = true;
                      last_label = spots[i];
                      continue;
                    }
                    edge.id1 = my_label;
                    edge.id2 = spots[i];
                    edge.pred = pred;
                    edge.boundary = false;
                    sorted_edges[sorted_edges_len++] = edge;
                    last_label = spots[i];
                  }
                  // add a self edge.
                  edge.id1 = my_label;
                  edge.id2 = my_label;
                  edge.pred = pred;
                  edge.boundary = (boundary || !my_label);
                  sorted_edges[sorted_edges_len++] = edge;
              }
          }
      }

    }
 
    std::sort (sorted_edges, sorted_edges + sorted_edges_len, tfk_compare);
    rag_add_edge_batch(rag, sorted_edges, sorted_edges_len, feature_man);
    free(sorted_edges);
}

void Stack::build_rag_loop_leek(Rag<Label>* rag, FeatureMgr* feature_man, 
				int x_start, int x_end, int y_start, int y_end,
				int z_start, int z_end, bool use_mito_prob)
{
    if (feature_man->get_num_channels() > 7) {
	throw ErrMsg("Cannot run with more than 7 channels");
    }
    
    unsigned int maxx = width-1; 
    unsigned int maxy = height-1; 
    unsigned int maxz = depth-1; 

    std::vector<leek_edge> sorted_edges;

    unsigned int * label_vol_array = watershed;
    unsigned int plane_size = width * height;
    uint64_t label_vol_stride[3] = {1,width,plane_size}; 
    leek_edge edge;

    if (x_start < 0) x_start = 0;
    if (y_start < 0) y_start = 0;
    if (z_start < 0) z_start = 0;

    bool avoid_conditionals = true;
    if (z_start - 1 < 0 || y_start -1 < 0 || x_start - 1 < 0 || z_end+1>maxz || y_end+1 > maxy || x_end+1 > maxx) {
      avoid_conditionals = false;
    }


    unsigned int spots[7];

    if (avoid_conditionals) {
      unsigned int* spots_ptr[7];
      int64_t curr_spot = x_start + y_start*width + z_start*plane_size;
      spots_ptr[0] = watershed + curr_spot;
      spots_ptr[1] = watershed + curr_spot-1;
      spots_ptr[2] = watershed + curr_spot+1;
      spots_ptr[3] = watershed + curr_spot-width;
      spots_ptr[4] = watershed + curr_spot+width;
      spots_ptr[5] = watershed + curr_spot-plane_size;
      spots_ptr[6] = watershed + curr_spot+plane_size;

      for (unsigned int z = z_start; z < z_end; z++) {
	  int z_spot = z * plane_size;
	  spots_ptr[5] += plane_size;
	  spots_ptr[6] += plane_size;
	  for (unsigned int y = y_start; y < y_end; y++) {
	      int y_spot = y * width;
	      spots_ptr[3] += width;
	      spots_ptr[4] += width;
	      for (unsigned int x = x_start; x < x_end; x++) {
		  curr_spot = x + y_spot + z_spot;
		  spots_ptr[0] += 1;
		  spots_ptr[1] += 1;
		  spots_ptr[2] += 1;
		  for (int i = 0; i < 7; i++) {
		    spots[i] = *(spots_ptr[i]);
		  }


		  unsigned char pred[7];
		  for (int idx=0; idx<feature_man->get_num_channels(); idx++) {
		      pred[idx] = prediction_array[idx][curr_spot];
		  }

		  int my_label = spots[0];

		  // sort.
		  for (int j = 0; j < 7; j++) {
		     int insert = spots[j];
		     int hole = j;
		     while (hole > 0 && insert > spots[hole-1]) {
		       spots[hole] = spots[hole-1];
		       hole--;
		     }
		     spots[hole] = insert;
		  }
	  
		  int last_label = my_label;
		  bool boundary = false;

		  for (int i = 0; i < 7; i++) {
		    if (spots[i] == last_label || spots[i] == my_label) continue;
		    if (!spots[i]) {
		      boundary = true;
		      last_label = spots[i];
		      continue;
		    }
		    edge.id1 = my_label;
		    edge.id2 = spots[i];
		    for (size_t idx=0; idx<feature_man->get_num_channels(); idx++) {
			edge.pred[idx] = pred[idx];
		    }
		    edge.boundary = false;
		    sorted_edges.push_back(edge);
		    last_label = spots[i];
		  }
		  // add a self edge.
		  edge.id1 = my_label;
		  edge.id2 = my_label;
		  for (size_t idx=0; idx<feature_man->get_num_channels(); idx++) {
		    edge.pred[idx] = pred[idx];
		  }

		  edge.boundary = (boundary || !my_label);
		  sorted_edges.push_back(edge);
	      }
	  }
      }
    } else {
      for (unsigned int z = z_start; z < z_end; z++) {
	  int z_spot = z * plane_size;
	  //if (z < 0 || z > maxz) continue;
	  for (unsigned int y = y_start; y < y_end; y++) {
	      //if (y < 0 || y > maxy) continue;
	      int y_spot = y * width;
	      for (unsigned int x = x_start; x < x_end; x++) {
		  //if (x < 0 || x > maxx) continue;
		  int64_t curr_spot = x + y_spot + z_spot;
		  spots[0] = watershed[curr_spot];
		  spots[1] = x >= 0 ? watershed[curr_spot-1] : 0;
		  spots[2] = x+1 <= maxx ? watershed[curr_spot+1] : 0;
		  spots[3] = y-1 >=0 ? watershed[curr_spot-width] : 0;
		  spots[4] = y+1 <= maxy ? watershed[curr_spot+width] : 0;
		  spots[5] = z-1 >= 0 ? watershed[curr_spot-plane_size] : 0;
		  spots[6] = z+1 <= maxz ? watershed[curr_spot+plane_size] : 0;

		  unsigned char pred[7];
		  for (int idx=0; idx<feature_man->get_num_channels(); idx++) {
		      pred[idx] = prediction_array[idx][curr_spot];
		  }

		  int my_label = spots[0];

		  // sort.
		  for (int j = 0; j < 7; j++) {
		     int insert = spots[j];
		     int hole = j;
		     while (hole > 0 && insert > spots[hole-1]) {
		       spots[hole] = spots[hole-1];
		       hole--;
		     }
		     spots[hole] = insert;
		  }
	  
		  int last_label = my_label;
		  bool boundary = false;

		  for (int i = 0; i < 7; i++) {
		    if (spots[i] == last_label || spots[i] == my_label) continue;
		    if (!spots[i]) {
		      boundary = true;
		      last_label = spots[i];
		      continue;
		    }
		    edge.id1 = my_label;
		    edge.id2 = spots[i];
		    for (size_t idx=0; idx<feature_man->get_num_channels(); idx++) {
			edge.pred[idx] = pred[idx];
		    }
		    edge.boundary = false;
		    sorted_edges.push_back(edge);
		    last_label = spots[i];
		  }
		  // add a self edge.
		  edge.id1 = my_label;
		  edge.id2 = my_label;
  		  for (size_t idx=0; idx<feature_man->get_num_channels(); idx++) {
		    edge.pred[idx] = pred[idx];
		  }
		  edge.boundary = (boundary || !my_label);
		  sorted_edges.push_back(edge);
	      }
	  }
      }

    }
 
    std::sort (sorted_edges.begin(), sorted_edges.end(), leek_compare);
    rag_add_edge_batch(rag, &sorted_edges.begin()[0], sorted_edges.size(), feature_man);
}

inline void move_edge_feature (FeatureMgr* fm1, FeatureMgr* fm2, RagEdge<Label>* edge1, RagEdge<Label>* edge2) {
    // cout << "Move Edge" << endl;
    edge1->set_size(edge2->get_size()); 
    EdgeCaches &ec1 = fm1->get_edge_cache();
    EdgeCaches &ec2 = fm2->get_edge_cache();
    EdgeCaches::iterator edge_feat2 = ec2.find(edge2);
    assert(edge_feat2 != ec2.end());
    ec1[edge1] = ec2[edge2];
    ec2[edge2] = std::vector<void *>();
    // fm1->add_val(0.0, edge1);
    // merge_edge_features(fm1, fm2, edge1, edge2);
}

void merge_edge_features (FeatureMgr* fm1, FeatureMgr* fm2, RagEdge<Label>* edge1, RagEdge<Label>* edge2) {
    // cout << "Merge Edge" << endl;
    edge1->incr_size(edge2->get_size());
    EdgeCaches &ec1 = fm1->get_edge_cache();
    EdgeCaches &ec2 = fm2->get_edge_cache();
    unsigned int pos = 0;
    unsigned int num_chan = fm1->get_num_channels();
    for (int i = 0; i < num_chan; ++i) {
        vector<FeatureCompute*>& features = fm1->get_channel_features()[i];
        for (int j = 0; j < features.size(); ++j) {
            if (ec1[edge1][pos] && ec2[edge2][pos]) {
                features[j]->merge_cache(ec1[edge1][pos], ec2[edge2][pos]);
            }
            ++pos;
        }
    }    
}

void move_node_feature (FeatureMgr* fm1, FeatureMgr* fm2, RagNode<Label>* node1, RagNode<Label>* node2) {
    // cout << "Move Node" << endl;
    node1->set_size(node2->get_size());
    node1->set_border_size(node2->get_border_size());
    NodeCaches &nc1 = fm1->get_node_cache();
    NodeCaches &nc2 = fm2->get_node_cache();
    NodeCaches::iterator node_feat2 = nc2.find(node2);
    // assert(node_feat2 != nc2.end());
    if (node_feat2 == nc2.end())
        return;
    nc1[node1] = nc2[node2];
    nc2[node2] = std::vector<void *>();
    // fm1->add_val(0.0, node1);
    // merge_node_features(fm1, fm2, node1, node2);
}

void merge_node_features (FeatureMgr* fm1, FeatureMgr* fm2, RagNode<Label>* node1, RagNode<Label>* node2) {
    // cout << "Merge Node" << endl;
    NodeCaches &nc1 = fm1->get_node_cache();
    NodeCaches &nc2 = fm2->get_node_cache();
    if (nc2.find(node2) == nc2.end()) {
        return;
    }
    if (nc1.find(node1) == nc1.end()) {
        move_node_feature(fm1, fm2, node1, node2);
        return;
    }
    node1->incr_size(node2->get_size());
    node1->incr_border_size(node2->get_border_size());
    unsigned int pos = 0;
    unsigned int num_chan = fm1->get_num_channels();
    // cout << "Node 2: " << node2->get_node_id() << endl;
    for (int i = 0; i < num_chan; ++i) {
        vector<FeatureCompute*>& features = fm1->get_channel_features()[i];
        for (int j = 0; j < features.size(); ++j) {
            if (nc1[node1][pos] && nc2[node2][pos]) {
                features[j]->merge_cache(nc1[node1][pos], nc2[node2][pos]);
            }
            ++pos;
        }
    }
}



void Stack::merge_rags (Rag<Label>* rag1, Rag<Label>* rag2, FeatureMgr* fm1, FeatureMgr* fm2) {
    set<Label> processed;
    for (Rag<Label>::nodes_iterator it1 = rag2->nodes_begin(); it1 != rag2->nodes_end(); ++it1) {
        assert(processed.find((*it1)->get_node_id()) == processed.end());
        RagNode<Label> * node1 = rag1->find_rag_node((*it1)->get_node_id());
        if (!node1) {
            // if not, insert the node and its incident edges
            RagNode<Label>* new_node = rag1->insert_rag_node((*it1)->get_node_id());
            move_node_feature(fm1, fm2, new_node, *it1);
            // assert(fm1->get_node_cache().find(new_node) != fm1->get_node_cache().end());
            for (RagNode<Label>::edge_iterator it2 = (*it1)->edge_begin(); it2 != (*it1)->edge_end(); ++it2) {
                RagNode<Label>* terminal_node = (*it2)->get_other_node(*it1);
                node1 = rag1->find_rag_node(terminal_node->get_node_id());
                if (node1) {
                    // add edge and update node
                    RagEdge<Label>* new_edge = rag1->insert_rag_edge(node1, new_node);
                    move_edge_feature(fm1, fm2, new_edge, *it2);
                    // assert(fm1->get_edge_cache().find(new_edge) != fm1->get_edge_cache().end());
                }
            }
        } else {
            // merge size. go thru neighbors. if not processed and in rag1, update edge.
            assert(processed.find(node1->get_node_id()) == processed.end());
            merge_node_features(fm1, fm2, node1, *it1);
            for (RagNode<Label>::edge_iterator it2 = (*it1)->edge_begin(); it2 != (*it1)->edge_end(); ++it2) {
                RagNode<Label>* terminal_node = (*it2)->get_other_node(*it1);
                RagNode<Label>* node2 = rag1->find_rag_node(terminal_node->get_node_id()); 
                if (node2) {
                    if (processed.find(node2->get_node_id()) == processed.end()) {
                        RagEdge<Label>* new_edge = rag1->find_rag_edge(node1->get_node_id(), node2->get_node_id());
                        if (new_edge) {
                            // if edge is already there
                            merge_edge_features(fm1, fm2, new_edge, *it2);
                        } else {
                            new_edge = rag1->insert_rag_edge(node1, node2);
                            move_edge_feature(fm1, fm2, new_edge, *it2);
                        }
                    }
                }
            }
        }
        processed.insert((*it1)->get_node_id());
    }
}


void Stack::build_rag_recurse (Rag<Label>* rag1, FeatureMgr* fm1,  
                                                    int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob) {
    int x_size = x_end - x_start;
    int y_size = y_end - y_start;
    int z_size = z_end - z_start;

    Rag<Label>* rag2 = new Rag<Label>();
    FeatureMgr* fm2 = new FeatureMgr();
    fm2->copy_channel_features(fm1);
    //FeatureMgr* fm2 = new FeatureMgr(prediction_array.size());
    //_set_basic_features(fm2); 

    bool stop_recurse = false;

    if (x_size >= y_size && x_size >= z_size) {
        if (x_size > BLOCK_SIZE_LIMIT) {
            cilk_spawn build_rag_recurse(rag1, fm1,  x_start, x_start + x_size/2, y_start, y_end, z_start, z_end, use_mito_prob);
            build_rag_recurse(rag2, fm2, x_start + x_size/2, x_end, y_start, y_end, z_start, z_end, use_mito_prob);
            cilk_sync;
        } else {
            stop_recurse = true;
        }
    } else if (y_size >= x_size && y_size >= z_size) {
        if (y_size > BLOCK_SIZE_LIMIT) {
            cilk_spawn build_rag_recurse(rag1, fm1, x_start, x_end, y_start, y_start + y_size/2, z_start, z_end, use_mito_prob);
            build_rag_recurse(rag2, fm2, x_start, x_end, y_start + y_size/2, y_end, z_start, z_end, use_mito_prob);
            cilk_sync;
        } else {
            stop_recurse = true;
        }
    } else {
        if (z_size > BLOCK_SIZE_LIMIT_Z) {
            cilk_spawn build_rag_recurse(rag1, fm1, x_start, x_end, y_start, y_end, z_start, z_start + z_size/2, use_mito_prob);
            build_rag_recurse(rag2, fm2, x_start, x_end, y_start, y_end, z_start + z_size/2, z_end, use_mito_prob);
            cilk_sync;
        } else {
            stop_recurse = true;
        }
    }

    if (stop_recurse) {
        build_rag_loop(rag1, fm1, x_start, x_end, y_start, y_end, z_start, z_end, use_mito_prob);
    } else {
        merge_rags(rag1, rag2, fm1, fm2);
        //if (use_mito_prob)
        //    merge_mito_probs (mito_prob1, mito_prob2);
    }
}




void Stack::build_rag()
{
    if (feature_mgr && (feature_mgr->get_num_features() == 0)) {
        feature_mgr->add_median_feature();
        median_mode = true;
        printf("In median mode wtf?\n");
    } 
    
    unsigned int plane_size = width * height;
    std::vector<double> predictions(prediction_array.size(), 0.0);
    std::tr1::unordered_set<Label> labels;
    unsigned long long curr_spot;
    
    unsigned int wcount=0;
    int x_start = 1;
    int x_end = width-1;
    int y_start = 1;
    int y_end = height-1;
    int z_start = 1;
    int z_end = depth-1;
    bool use_mito_prob = false; 
    if (__cilkrts_get_nworkers() > 1) {
	build_rag_recurse (rag, feature_mgr, x_start, x_end, y_start, y_end, z_start, z_end, false);
    } else {
	build_rag_loop (rag, feature_mgr, x_start, x_end, y_start, y_end, z_start, z_end, false);
    }

    if (merge_mito)
	printf("Deciding mito sps with thd= %.3lf ...", mito_thd);
    watershed_to_body[0] = 0;
    for (Rag<Label>::nodes_iterator iter = rag->nodes_begin(); iter != rag->nodes_end(); ++iter) {
        Label id = (*iter)->get_node_id();
        watershed_to_body[id] = id;
	if (merge_mito)
	    (*iter)->get_type_decider()->set_type(mito_thd);
    }

    
     //unsigned int ecount = 0;
     //for (Rag<Label>::edges_iterator iter = rag->edges_begin(); iter != rag->edges_end(); ++iter) {
     //   (*iter)->set_edge_id(ecount);
     //   ecount++;     
     //}
}




int Stack::remove_inclusions()
{
    int num_removed = 0;

    visited.clear();
    node_depth.clear();
    low_count.clear();
    prev_id.clear();
    biconnected_components.clear();
    stack.clear();

    RagNode<Label>* rag_node = 0;
    for (Rag<Label>::nodes_iterator iter = rag->nodes_begin(); iter != rag->nodes_end(); ++iter) {
        bool border = (*iter)->is_border();

        if (border) {
            rag_node = *iter;
            break;
        }
    }
    assert(rag_node);

    DFSStack temp;
    temp.previous = 0;
    temp.rag_node = rag_node;
    temp.count = 1;
    temp.start_pos = 0;

    std::vector<DFSStack> dfs_stack;
    dfs_stack.push_back(temp);
    biconnected_dfs(dfs_stack);
    stack.clear();
    std::tr1::unordered_map<Label, Label> body_to_body; 
    std::tr1::unordered_map<Label, std::vector<Label> > merge_history2; 

    // merge nodes in biconnected_components (ignore components with '0' node)
    for (int i = 0; i < biconnected_components.size(); ++i) {
        bool found_zero = false;
        std::tr1::unordered_set<Label> merge_nodes;
        for (int j = 0; j < biconnected_components[i].size()-1; ++j) {
            Label region1 = biconnected_components[i][j].region1;     
            Label region2 = biconnected_components[i][j].region2;
            if (!region1 || !region2) {
                found_zero = true;
            }
            
            if (body_to_body.find(region1) != body_to_body.end()) {
                region1 = body_to_body[region1];
            }
            if (body_to_body.find(region2) != body_to_body.end()) {
                region2 = body_to_body[region2];
            }

            assert(region1 != region2);

            merge_nodes.insert(region1);
            merge_nodes.insert(region2);
        }
        if (!found_zero) {
            Label articulation_region = biconnected_components[i][biconnected_components[i].size()-1].region1;
            unsigned long long total_size = 0;
            RagNode<Label>* articulation_node = rag->find_rag_node(articulation_region);
            for (std::tr1::unordered_set<Label>::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter) {
                Label region1 = *iter;
                RagNode<Label>* rag_node = rag->find_rag_node(region1);
                total_size += rag_node->get_size();

            }
          
            bool found_preserve = false; 
            for (std::tr1::unordered_set<Label>::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter) {
                Label region2 = *iter;
                if (body_to_body.find(region2) != body_to_body.end()) {
                    region2 = body_to_body[region2];
                }

                RagNode<Label>* rag_node = rag->find_rag_node(region2);
                assert(rag_node);
                if (articulation_node != rag_node) {
                    for (RagNode<Label>::edge_iterator edge_iter = rag_node->edge_begin();
                        edge_iter != rag_node->edge_end(); ++edge_iter) {
                        if ((*edge_iter)->is_preserve()) {
                            found_preserve = true;
                            break;
                        }
                    } 
                }

                if (found_preserve) {
                    break;
                }
            }

            if (found_preserve) {
                continue;
            }
            
            for (std::tr1::unordered_set<Label>::iterator iter = merge_nodes.begin(); iter != merge_nodes.end(); ++iter) {
                Label region2 = *iter;
                if (body_to_body.find(region2) != body_to_body.end()) {
                    region2 = body_to_body[region2];
                }

                RagNode<Label>* rag_node = rag->find_rag_node(region2);
                assert(rag_node);
                if (articulation_node != rag_node) {
                    feature_mgr->merge_features(articulation_node, rag_node); 
                    rag->remove_rag_node(rag_node); 

                    watershed_to_body[region2] = articulation_region;
                    for (std::vector<Label>::iterator iter2 = merge_history[region2].begin(); iter2 != merge_history[region2].end(); ++iter2) {
                        watershed_to_body[*iter2] = articulation_region;
                    }

                    body_to_body[region2] = articulation_region;
                    for (std::vector<Label>::iterator iter2 = merge_history2[region2].begin(); iter2 != merge_history2[region2].end(); ++iter2) {
                        body_to_body[*iter2] = articulation_region;
                    }

                    merge_history[articulation_region].push_back(region2);
                    merge_history[articulation_region].insert(merge_history[articulation_region].end(), merge_history[region2].begin(), merge_history[region2].end());
                    merge_history.erase(region2);
                   
                    merge_history2[articulation_region].push_back(region2); 
                    merge_history2[articulation_region].insert(merge_history2[articulation_region].end(), merge_history2[region2].begin(), merge_history2[region2].end());
                    merge_history2.erase(region2);
                }
            }
        }
    }

    
    visited.clear();
    node_depth.clear();
    low_count.clear();
    prev_id.clear();
    biconnected_components.clear();
    stack.clear();
    return num_removed;
}

// return low point (0 is default edge)
void Stack::biconnected_dfs(std::vector<DFSStack>& dfs_stack)
{
    while (!dfs_stack.empty()) {
        DFSStack entry = dfs_stack.back();
        RagNode<Label>* rag_node = entry.rag_node;
        Label previous = entry.previous;
        int count = entry.count;
        dfs_stack.pop_back();

        if (visited.find(rag_node->get_node_id()) == visited.end()) {
            visited.insert(rag_node->get_node_id());
            node_depth[rag_node->get_node_id()] = count;
            low_count[rag_node->get_node_id()] = count;
            prev_id[rag_node->get_node_id()] = previous;
        }

        bool skip = false;
        int curr_pos = 0;
        for (RagNode<Label>::node_iterator iter = rag_node->node_begin(); iter != rag_node->node_end(); ++iter) {
            RagEdge<Label>* rag_edge = rag->find_rag_edge(rag_node, *iter);
            if (rag_edge->is_false_edge()) {
                continue;
            }
            
            if (curr_pos < entry.start_pos) {
                ++curr_pos;
                continue;
            }
            if (prev_id[(*iter)->get_node_id()] == rag_node->get_node_id()) {
                OrderedPair<Label> current_edge(rag_node->get_node_id(), (*iter)->get_node_id());
                int temp_low = low_count[(*iter)->get_node_id()];
                low_count[rag_node->get_node_id()] = std::min(low_count[rag_node->get_node_id()], temp_low);

                if (temp_low >= count) {
                    OrderedPair<Label> popped_edge;
                    biconnected_components.push_back(std::vector<OrderedPair<Label> >());
                    do {
                        popped_edge = stack.back();
                        stack.pop_back();
                        biconnected_components[biconnected_components.size()-1].push_back(popped_edge);
                    } while (!(popped_edge == current_edge));
                    OrderedPair<Label> articulation_pair(rag_node->get_node_id(), rag_node->get_node_id());
                    biconnected_components[biconnected_components.size()-1].push_back(articulation_pair);
                } 
            } else if (visited.find((*iter)->get_node_id()) == visited.end()) {
                OrderedPair<Label> current_edge(rag_node->get_node_id(), (*iter)->get_node_id());
                stack.push_back(current_edge);

                DFSStack temp;
                temp.previous = rag_node->get_node_id();
                temp.rag_node = (*iter);
                temp.count = count+1;
                temp.start_pos = 0;
                entry.start_pos = curr_pos;
                dfs_stack.push_back(entry);
                dfs_stack.push_back(temp);
                skip = true;
                break;
            } else if ((*iter)->get_node_id() != previous) {
                low_count[rag_node->get_node_id()] = std::min(low_count[rag_node->get_node_id()], node_depth[(*iter)->get_node_id()]);
                if (count > node_depth[(*iter)->get_node_id()]) {
                    stack.push_back(OrderedPair<Label>(rag_node->get_node_id(), (*iter)->get_node_id()));
                }
            }
            ++curr_pos;
        }

        if (skip) {
            continue;
        }

        bool border = rag_node->is_border();
        if (previous && border) {
            low_count[rag_node->get_node_id()] = 0;
            stack.push_back(OrderedPair<Label>(0, rag_node->get_node_id()));
        }
    }
}

Label* Stack::get_label_volume(){

        Label * temp_labels = new Label[(width-(2*padding))*(height-(2*padding))*(depth-(2*padding))];
        Label * temp_labels_iter = temp_labels;
        unsigned int plane_size = width * height;

        for (int z = padding; z < (depth - padding); ++z) {
            int z_spot = z * plane_size;
            for (int y = padding; y < (height - padding); ++y) {
                int y_spot = y * width;
                for (int x = padding; x < (width - padding); ++x) {
                    unsigned long long curr_spot = x+y_spot+z_spot;
                    Label sv_id = watershed[curr_spot];
                    Label body_id;
                    if (!watershed_to_body.empty()) {
                        body_id = watershed_to_body[sv_id];
                    } else {
                        body_id = sv_id;
                    }
                    *temp_labels_iter = body_id;
                    ++temp_labels_iter;
                }
            }
        }
    
        return temp_labels;
}






/*

void Stack::determine_edge_locations(){
        best_edge_z.clear();
        best_edge_loc.clear();

        unsigned int plane_size = width * height;

        for (unsigned int z = 1; z < (depth-1); ++z) {
            int z_spot = z * plane_size;

            EdgeCount curr_edge_z;
            EdgeLoc curr_edge_loc;
             

            for (unsigned int y = 1; y < (height-1); ++y) {
                int y_spot = y * width;
                for (unsigned int x = 1; x < (width-1); ++x) {
                    unsigned long long curr_spot = x + y_spot + z_spot;
                    unsigned int spot0 = watershed_to_body[(watershed[curr_spot])];
                    unsigned int spot1 = watershed_to_body[(watershed[curr_spot-1])];
                    unsigned int spot2 = watershed_to_body[(watershed[curr_spot+1])];
                    unsigned int spot3 = watershed_to_body[(watershed[curr_spot-width])];
                    unsigned int spot4 = watershed_to_body[(watershed[curr_spot+width])];
                    unsigned int spot5 = watershed_to_body[(watershed[curr_spot-plane_size])];
                    unsigned int spot6 = watershed_to_body[(watershed[curr_spot+plane_size])];

                    if (spot1 && (spot0 != spot1)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot1);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                    if (spot2 && (spot0 != spot2)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot2);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                    if (spot3 && (spot0 != spot3)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot3);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                    if (spot4 && (spot0 != spot4)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot4);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                    if (spot5 && (spot0 != spot5)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot5);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                    if (spot6 && (spot0 != spot6)) {
                        RagEdge<Label>* edge = rag->find_rag_edge(spot0, spot6);
                        curr_edge_z[edge] += 1;  
                        curr_edge_loc[edge] = Location(x,y,z);  
                    }
                }
            }
        
            for (EdgeCount::iterator iter = curr_edge_z.begin(); iter != curr_edge_z.end(); ++iter) {
                if (iter->second > best_edge_z[iter->first]) {
                    best_edge_z[iter->first] = iter->second;
                    best_edge_loc[iter->first] = curr_edge_loc[iter->first];
                }
            }
        } 
        
        
        return;
}

bool Stack::add_edge_constraint(boost::python::tuple loc1, boost::python::tuple loc2) {
        Label body1 = get_body_id(boost::python::extract<int>(loc1[0]), boost::python::extract<int>(loc1[1]),
                                boost::python::extract<int>(loc1[2]));    
        Label body2 = get_body_id(boost::python::extract<int>(loc2[0]), boost::python::extract<int>(loc2[1]),
                                boost::python::extract<int>(loc2[2]));    

        if (body1 == body2) {
            return false;
        }

        RagEdge<Label>* edge = rag->find_rag_edge(body1, body2);

        if (!edge) {
            RagNode<Label>* node1 = rag->find_rag_node(body1);
            RagNode<Label>* node2 = rag->find_rag_node(body2);
            edge = rag->insert_rag_edge(node1, node2);
            edge->set_false_edge(true);
        }
        edge->set_preserve(true);

        return true;
}
*/
Label Stack::get_body_id(unsigned int x, unsigned int y, unsigned int z) {
        x += padding;
        //y += padding;
        y = height - y - 1 - padding;
        z += padding;
        unsigned int plane_size = width * height;
        unsigned long long curr_spot = x + y * width + z * plane_size; 
        Label body_id = watershed[curr_spot];
        if (!watershed_to_body.empty()) {
            body_id = watershed_to_body[body_id];
        }
        return body_id;
}

/*boost::python::tuple Stack::get_edge_loc(RagEdge<Label>* edge) {
        if (best_edge_loc.find(edge) == best_edge_loc.end()) {
            throw ErrMsg("Edge location was not loaded!");
        }
        Location loc = best_edge_loc[edge];
        unsigned int x = boost::get<0>(loc) - padding;
        unsigned int y = boost::get<1>(loc) - padding; //height - boost::get<1>(loc) - 1 - padding;
        unsigned int z = boost::get<2>(loc) - padding;

        return boost::python::make_tuple(x, y, z); 
}

void Stack::get_edge_loc(RagEdge<Label>* edge, Label&x, Label& y, Label& z) {
        if (best_edge_loc.find(edge) == best_edge_loc.end()) {
            throw ErrMsg("Edge location was not loaded!");
        }
        Location loc = best_edge_loc[edge];
        x = boost::get<0>(loc) - padding;
        y = boost::get<1>(loc) - padding; //height - boost::get<1>(loc) - 1 - padding;
        z = boost::get<2>(loc) - padding;

}

void Stack::write_graph(string output_path){
  
    unsigned found = output_path.find_last_of(".");
    string dirname = output_path.substr(0,found);
    string filename = dirname + string("_edge_list.txt");
    
    FILE *fp = fopen(filename.c_str(), "wt");
    if(!fp){
        printf("%s could not be opened \n", filename.c_str());
	return;
    }
    fprintf(fp,"Nodes:\n");
    for (Rag<Label>::nodes_iterator iter = rag->nodes_begin(); iter != rag->nodes_end(); ++iter) {
	fprintf(fp,"%u %lu %u\n", (*iter)->get_node_id(), (*iter)->get_size(), 
		  (is_orphan(*iter)? 1 : 0));
    }
    fprintf(fp,"\nEdges:\n");
    
    Label x, y, z;
    for (Rag<Label>::edges_iterator iter = rag->edges_begin(); iter != rag->edges_end(); ++iter) {
        get_edge_loc(*iter, x, y, z); 
        fprintf(fp,"%u %u %lu %lu %u %u %lf %u %u %u\n", (*iter)->get_node1()->get_node_id(), (*iter)->get_node2()->get_node_id(),
	      (*iter)->get_node1()->get_size(), (*iter)->get_node2()->get_size(),
	      (*iter)->is_preserve()? 1:0, (*iter)->is_false_edge()? 1: 0,
	      (*iter)->get_weight(), x, y, z);
    }
    fclose(fp);
    
}

*/



int Stack::decide_edge_label(RagNode<Label>* rag_node1, RagNode<Label>* rag_node2){
    int edge_label = 0;
    Label node1 = rag_node1->get_node_id();	
    Label node2 = rag_node2->get_node_id();	

    Label node1gt = assignment.find(node1)->second;
    Label node2gt = assignment.find(node2)->second;

    if (!node1gt){
	printf("no grountruth node %d\n",node1);
	return 0;
    }  	
    if (!node2gt){
	printf("no grountruth node %d\n",node2);
	return 0;
    }  	

    edge_label = (node1gt == node2gt)? -1 : 1;

    if(merge_mito){
	if ( rag_node1->get_node_type() == 2 && rag_node2->get_node_type() == 2 ){
	    return 0;
	}
	else if ( rag_node1->get_node_type() == 2 || rag_node2->get_node_type() == 2 )
	    edge_label = 1;				
    }

    return edge_label;	
}





void Stack::set_basic_features(){

    double hist_percentiles[]={0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> percentiles(hist_percentiles, hist_percentiles+sizeof(hist_percentiles)/sizeof(double));		

    // ** for Toufiq's version ** feature_mgr->add_inclusiveness_feature(true);  	
    feature_mgr->add_moment_feature(4,true);	
    feature_mgr->add_hist_feature(25,percentiles,false); 	

    //cout << "Number of features, channels:" << feature_mgr->get_num_features()<< ","<<feature_mgr->get_num_channels()<<"\n";	

}









