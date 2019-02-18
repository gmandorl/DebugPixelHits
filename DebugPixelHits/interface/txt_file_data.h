#ifndef txt_file_data_h
#define txt_file_data_h

class txt_file_data {
      public:
          
          int index;
          
          int run;
          int lumi;
          int bx;
          int instLumi;
          int npv;
          int detid;
          
          float track_pt;
          float pred_hit_pos_y;
          float hit_pos_y;
          float cluster_center_y;
          float abs_hit_pos_y;
          float hit_pos_x;
          float cluster_center_x;
          float abs_hit_pos_x;
          float abs_pred_hit_pos_y;
          
          void set_index(int); 
          
          void set_run(int);
          void set_lumi(int);
          void set_bx(int);
          void set_instLumi(int);
          void set_npv(int);
          void set_detid(int);
          
          
          void set_track_pt(float);
          void set_pred_hit_pos_y(float);
          void set_hit_pos_y(float);
          void set_cluster_center_y(float);
          void set_abs_hit_pos_y(float);
          void set_hit_pos_x(float);
          void set_cluster_center_x(float);
          void set_abs_hit_pos_x(float);
          void set_abs_pred_hit_pos_y(float);
          
          
          
};


inline void txt_file_data::set_run(int a){run=a;};
inline void txt_file_data::set_lumi(int a){lumi=a;};
inline void txt_file_data::set_bx(int a){bx=a;};
inline void txt_file_data::set_instLumi(int a){instLumi=a;};
inline void txt_file_data::set_npv(int a){npv=a;};
inline void txt_file_data::set_detid(int a){detid=a;}; 
inline void txt_file_data::set_index(int a){index=a;};
inline void txt_file_data::set_track_pt(float a){track_pt=a;};
inline void txt_file_data::set_pred_hit_pos_y(float a){pred_hit_pos_y=a;};
inline void txt_file_data::set_hit_pos_y(float a){hit_pos_y=a;};
inline void txt_file_data::set_cluster_center_y(float a){cluster_center_y=a;};
inline void txt_file_data::set_abs_hit_pos_y(float a){abs_hit_pos_y=a;};
inline void txt_file_data::set_hit_pos_x(float a){hit_pos_x=a;};
inline void txt_file_data::set_cluster_center_x(float a){cluster_center_x=a;};
inline void txt_file_data::set_abs_hit_pos_x(float a){abs_hit_pos_x=a;};
inline void txt_file_data::set_abs_pred_hit_pos_y(float a){abs_pred_hit_pos_y=a;};

#endif
