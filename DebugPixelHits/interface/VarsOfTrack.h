#ifndef VarsOfTrack_h
#define VarsOfTrack_h

class VarsOfTrack {
      public:
          
                
        bool found;
        int detid, detIsActive, roc;
        float chargeCorr;
        float track_global_phi,track_global_z;
        float track_local_x,track_local_y;
        float track_local_Dx,track_local_Dy;
        float track_localPixel_x,track_localPixel_y;
        float hit_localPixel_x,hit_localPixel_y;
        float cluster_localPixel_x,cluster_localPixel_y;
        float brokenCluster_localPixel_x,brokenCluster_localPixel_y;
        float track_exp_sizeX, track_exp_sizeY, track_exp_charge;        
        float hit_global_phi, hit_global_z;
        float hit_local_x, hit_local_y;
        float cluster_center_x, cluster_center_y;
        float hit_firstpixel_x, hit_firstpixel_y;
        float hit_chi2, hit_charge;
        float hit_sizeX, hit_sizeY;
        float hitFound, hitOnTrack; 
        int hitInRandomWindow;
        float hitInRandomWindowDistance;
        int source_det,source_layer;        
        int maybeBadROC, trackHasHit, trackHasLostHit;
        float alpha, beta;        
        int cluster_charge_in_hits[147];
        int cluster_col_status[147];
        int cluster_chargeBroken_in_hits[147];
        int cluster_colBroken_status[147];
        
        bool column1_has_hit[416];
        bool column1_status[416];        
        bool column2_has_hit[416];
        bool column2_status[416];
        
        std::vector<int> All_hits_charge;
        std::vector<int> All_hits_Px;
        std::vector<int> All_hits_Py;
        
        
        bool broken_cl_FLAG, corr_FLAG;
        
        float dz_wCorrection, dz_compare;
        
        void set_found(bool);
        void set_detid(int);
        void set_detIsActive(int);
        void set_roc(int);
        void set_chargeCorr(float);
        void set_track_global_phi(float);
        void set_track_global_z(float);
        void set_track_local_x(float);
        void set_track_local_y(float);
        void set_track_local_Dx(float);
        void set_track_local_Dy(float);
        void set_track_localPixel_x(float);
        void set_track_localPixel_y(float);
        void set_hit_localPixel_x(float);
        void set_hit_localPixel_y(float);
        void set_cluster_localPixel_x(float);
        void set_cluster_localPixel_y(float);
        void set_brokenCluster_localPixel_x(float);
        void set_brokenCluster_localPixel_y(float);
        void set_track_exp_sizeX(float); 
        void set_track_exp_sizeY(float); 
        void set_track_exp_charge(float);        
        void set_hit_global_phi(float); 
        void set_hit_global_z(float);
        void set_hit_local_x(float); 
        void set_hit_local_y(float);
        void set_cluster_center_x(float); 
        void set_cluster_center_y(float);
        void set_hit_firstpixel_x(float);
        void set_hit_firstpixel_y(float);
        void set_hit_chi2(float);
        void set_hit_charge(float);
        void set_hit_sizeX(float);
        void set_hit_sizeY(float);
        void set_hitFound(float);
        void set_hitOnTrack(float); 
        void set_hitInRandomWindow(int);
        void set_hitInRandomWindowDistance(float);
        void set_source_det(float);
        void set_source_layer(float);        
        void set_maybeBadROC(int); 
        void set_trackHasHit(int); 
        void set_trackHasLostHit(int);         
        void set_alpha(float);
        void set_beta(float);        
        void set_cluster_charge_in_hits(int[147] );
        void set_cluster_col_status(int[147] );
        void set_cluster_chargeBroken_in_hits(int[147] );
        void set_cluster_colBroken_status(int[147] );
        
        void set_column1_has_hit(bool[416]);
        void set_column1_status(bool[416]);        
        void set_column2_has_hit(bool[416]);
        void set_column2_status(bool[416]);    
        
        void set_All_hits_charge(std::vector<int>);
        void set_All_hits_Px(std::vector<int>);
        void set_All_hits_Py(std::vector<int>);
        
        void set_broken_Cluster_Flag(bool);
        void set_Correction_Flag(bool);
        void set_dz_wCorrection(float);
        void set_dz_compare(float);
                
        int example;        
        void set_example(int);
        
           
          
};

inline void VarsOfTrack::set_found(bool a){found=a;}
inline void VarsOfTrack::set_detid(int a){detid=a;}
inline void VarsOfTrack::set_detIsActive(int a){detIsActive=a;}
inline void VarsOfTrack::set_roc(int a){roc=a;}
inline void VarsOfTrack::set_chargeCorr(float a){chargeCorr=a;}
inline void VarsOfTrack::set_track_global_phi(float a){track_global_phi=a;}
inline void VarsOfTrack::set_track_global_z(float a){track_global_z=a;}
inline void VarsOfTrack::set_track_local_x(float a){track_local_x=a;}
inline void VarsOfTrack::set_track_local_y(float a){track_local_y=a;}
inline void VarsOfTrack::set_track_local_Dy(float a){track_local_Dy=a;}
inline void VarsOfTrack::set_track_local_Dx(float a){track_local_Dx=a;}
inline void VarsOfTrack::set_track_localPixel_x(float a){track_localPixel_x=a;}
inline void VarsOfTrack::set_track_localPixel_y(float a){track_localPixel_y=a;}
inline void VarsOfTrack::set_hit_localPixel_x(float a){hit_localPixel_x=a;}
inline void VarsOfTrack::set_hit_localPixel_y(float a){hit_localPixel_y=a;}
inline void VarsOfTrack::set_cluster_localPixel_x(float a){cluster_localPixel_x=a;}
inline void VarsOfTrack::set_cluster_localPixel_y(float a){cluster_localPixel_y=a;}
inline void VarsOfTrack::set_brokenCluster_localPixel_x(float a){brokenCluster_localPixel_x=a;}
inline void VarsOfTrack::set_brokenCluster_localPixel_y(float a){brokenCluster_localPixel_y=a;}
inline void VarsOfTrack::set_track_exp_sizeX(float a){track_exp_sizeX=a;} 
inline void VarsOfTrack::set_track_exp_sizeY(float a){track_exp_sizeY=a;} 
inline void VarsOfTrack::set_track_exp_charge(float a){track_exp_charge=a;}        
inline void VarsOfTrack::set_hit_global_phi(float a){hit_global_phi=a;} 
inline void VarsOfTrack::set_hit_global_z(float a){hit_global_z=a;}
inline void VarsOfTrack::set_hit_local_x(float a){hit_local_x=a;} 
inline void VarsOfTrack::set_hit_local_y(float a){hit_local_y=a;}
inline void VarsOfTrack::set_cluster_center_x(float a){cluster_center_x=a;} 
inline void VarsOfTrack::set_cluster_center_y(float a){cluster_center_y=a;}
inline void VarsOfTrack::set_hit_firstpixel_x(float a){hit_firstpixel_x=a;}
inline void VarsOfTrack::set_hit_firstpixel_y(float a){hit_firstpixel_y=a;}
inline void VarsOfTrack::set_hit_chi2(float a){hit_chi2=a;}
inline void VarsOfTrack::set_hit_charge(float a){hit_charge=a;}
inline void VarsOfTrack::set_hit_sizeX(float a){hit_sizeX=a;}
inline void VarsOfTrack::set_hit_sizeY(float a){hit_sizeY=a;}
inline void VarsOfTrack::set_hitFound(float a){hitFound=a;}
inline void VarsOfTrack::set_hitOnTrack(float a){hitOnTrack=a;} 
inline void VarsOfTrack::set_hitInRandomWindow(int a){hitInRandomWindow=a;}
inline void VarsOfTrack::set_hitInRandomWindowDistance(float a){hitInRandomWindowDistance=a;}
inline void VarsOfTrack::set_source_det(float a){source_det=a;}
inline void VarsOfTrack::set_source_layer(float a){source_layer=a;}        
inline void VarsOfTrack::set_maybeBadROC(int a){maybeBadROC=a;} 
inline void VarsOfTrack::set_trackHasHit(int a){trackHasHit=a;} 
inline void VarsOfTrack::set_trackHasLostHit(int a){trackHasLostHit=a;}         
inline void VarsOfTrack::set_alpha(float a){alpha=a;}
inline void VarsOfTrack::set_beta(float a){beta=a;}
inline void VarsOfTrack::set_cluster_charge_in_hits(int a[147]){for(unsigned int i=0; i<147; i++){cluster_charge_in_hits[i]=a[i];}}
inline void VarsOfTrack::set_cluster_col_status(int a[147]){for(unsigned int i=0; i<147; i++){cluster_col_status[i]=a[i];}}
inline void VarsOfTrack::set_cluster_chargeBroken_in_hits(int a[147]){for(unsigned int i=0; i<147; i++){cluster_chargeBroken_in_hits[i]=a[i];}}
inline void VarsOfTrack::set_cluster_colBroken_status(int a[147]){for(unsigned int i=0; i<147; i++){cluster_colBroken_status[i]=a[i];}}

inline void VarsOfTrack::set_column1_has_hit(bool a[416]){for(unsigned int i=0; i<416; i++){column1_has_hit[i]=a[i];}};
inline void VarsOfTrack::set_column1_status(bool a[416]){for(unsigned int i=0; i<416; i++){column1_status[i]=a[i];}};        
inline void VarsOfTrack::set_column2_has_hit(bool a[416]){for(unsigned int i=0; i<416; i++){column2_has_hit[i]=a[i];}};
inline void VarsOfTrack::set_column2_status(bool a[416]){for(unsigned int i=0; i<416; i++){column2_status[i]=a[i];}};   

inline void VarsOfTrack::set_All_hits_charge(std::vector<int> a){for(unsigned int i=0; i<a.size(); i++){All_hits_charge.push_back(a[i]);}};
inline void VarsOfTrack::set_All_hits_Px(std::vector<int> a){for(unsigned int i=0; i<a.size(); i++){All_hits_Px.push_back(a[i]);}};
inline void VarsOfTrack::set_All_hits_Py(std::vector<int> a){for(unsigned int i=0; i<a.size(); i++){All_hits_Py.push_back(a[i]);}};

inline void VarsOfTrack::set_broken_Cluster_Flag(bool a){broken_cl_FLAG=a;}
inline void VarsOfTrack::set_Correction_Flag(bool a){corr_FLAG=a;}
inline void VarsOfTrack::set_example(int a){example=a;}

inline void VarsOfTrack::set_dz_wCorrection(float a){dz_wCorrection=a;}
inline void VarsOfTrack::set_dz_compare(float a){dz_compare=a;}
#endif
