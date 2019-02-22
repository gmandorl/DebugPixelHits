import numpy
seed = 7 
numpy.random.seed(seed)

from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn.metrics import roc_curve, auc


def evaluate(keras_model, npz_file, epoch):
    
    plt.clf()
    
    name_String="Eval_epoch"+str(epoch)+"_"
    
    print "\n"
    
    thefile=numpy.load(npz_file)
    sc = StandardScaler().fit(thefile["track"])
    print sc.mean_

    sc2 = StandardScaler().fit(thefile["track1"])
    print sc2.mean_
    print sc2.scale_
    
    len2=len(thefile["track"])

    range1=range(0,len2,2)
    range2=range(1,len2,2)
    
    ll=keras_model.predict([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]])
    
    
    #variabili del cluster
    ll2=thefile["hit_pos"][range2,1]
    delta=thefile["delta_TH"][range2,1]
    track_pos=thefile["track_pos"][range2,1]
    clucenter=thefile["local_cluster"][range2,1]
    cluBroken=thefile["local_brokencluster"][range2,1]
    broekn=thefile["Flag"][range2,0]
    
    print ll[0].shape, ll[1].shape, ll2.shape
    
    for i in range(100):
        print ll[0][i], "vs", ll2[i]
        #print thefile["truth1"][i,1]
        

    plt.hist(ll[0][:,0], bins=100, range=[-10,10], alpha=0.5)
    plt.hist(ll2, bins=100, range=[-10,10], alpha=0.5)
    plt.savefig(name_String+"plot.png")
    plt.clf()

    plt.hist(ll[0][:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_EstimateOK.png")
    plt.clf() 


    plt.hist(ll[0][:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.hist(delta, bins=100, range=[-2,2], alpha=0.5)
    plt.savefig(name_String+"plot_resolutionTrack1.png")
    plt.clf() 

    plt.hist(ll[0][:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(delta, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_resolutionTrack.png")
    plt.clf() 


    plt.hist(ll[1][:,0][thefile["Flag"][range2,0]==1], bins=100, range=[0,1], alpha=0.5, color="magenta", edgecolor='magenta')
    plt.hist(ll[1][:,0][thefile["Flag"][range2,0]==0], bins=100, range=[0,1], alpha=0.5, color="gold", edgecolor="gold")
    plt.savefig(name_String+"plot_TrueFalse.png")
    plt.clf() 

    #roc curve to be added--> take roc from traininig
    
    fpr, tpr, _ = roc_curve(thefile["Flag"][range2,0]==1, ll[1][:,0])    
    print "ROC AUC:  ",auc(fpr, tpr)
    
    plt.semilogy(tpr, fpr, label='ROC curve')       
    plt.semilogy([0, 1], [0, 1], 'k--')
    plt.xticks( [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] )
    plt.grid(True, which='both')
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('broken cluster efficiency')
    plt.ylabel('broken cluster mistag rate')
    plt.title('Roc curve, auc '+str(round(auc(fpr, tpr),2))) 
    #out_file.close()
    plt.savefig(name_String+"rocCurve_TrueFalse.png")
    plt.clf()
    

    plt.hist2d(ll2, ll[0][:,0],bins=100, norm=LogNorm())
    cb=plt.colorbar()
    plt.savefig(name_String+"plot_corrRainbowLog.png")
    cb.remove()
    plt.clf()


    plt.hist(ll[0][:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(ll[0][:,0]-clucenter, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(ll[0][:,0]-cluBroken, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_resolutionCluster.png")
    plt.clf()


    plt.hist(ll[0][:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken*(broekn==1)+clucenter*(broekn==0)-ll2, bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos.png")
    plt.clf()

    plt.hist(ll[0][broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyBroken.png")
    plt.clf()

    plt.hist(ll[0][broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyFull.png")
    plt.clf()

    plt.hist((ll[0][:,0]*(ll[1][:,0])+ll2*(1-ll[1][:,0]))-ll2, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken*(broekn==1)+clucenter*(broekn==0)-ll2, bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_Weighted.png")
    plt.clf()

    print ll[1][22,0]>0.7, "teu"

    for i in range(100):
        print (ll[0][:,0]*(ll[1][:,0]>0.7)+ll2*(ll[1][:,0]<=0.7))[i], ll[1][i,0], ll[0][i,0], ll2[i]

    plt.hist((ll[0][:,0]*(ll[1][:,0]>0.7)+ll2*(ll[1][:,0]<=0.7))-ll2, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken*(broekn==1)+clucenter*(broekn==0)-ll2, bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_Cut.png")
    plt.clf()


    plt.hist((ll[0][:,0]*(ll[1][:,0]>0.7)+ll2*(ll[1][:,0]<=0.7))[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyBrokenCUT.png")
    plt.clf()


    plt.hist((ll[0][:,0]*(ll[1][:,0])+clucenter*(1-ll[1][:,0]))-clucenter, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    #plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    #plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(ll[0][:,0]-(cluBroken*(broekn==1)+clucenter*(broekn==0)), bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_WeightedCluster.png")
    plt.clf()
        
        
def evaluateONLY_TrueFalse(keras_model, npz_file, epoch):
    
    plt.clf()
    
    name_String="Eval_epoch"+str(epoch)+"_"
    
    print "\n"
    
    thefile=numpy.load(npz_file)
    sc = StandardScaler().fit(thefile["track"])
    print sc.mean_

    sc2 = StandardScaler().fit(thefile["track1"])
    print sc2.mean_
    print sc2.scale_
    
    len2=len(thefile["track"])

    range1=range(0,len2,2)
    range2=range(1,len2,2)
    
    ll=keras_model.predict([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]])
    
    
    #variabili del cluster
    ll2=thefile["hit_pos"][range2,1]
    delta=thefile["delta_TH"][range2,1]
    track_pos=thefile["track_pos"][range2,1]
    clucenter=thefile["local_cluster"][range2,1]
    cluBroken=thefile["local_brokencluster"][range2,1]
    broekn=thefile["Flag"][range2,0]
    
    print ll.shape, ll2.shape
    
    
    plt.hist(ll[:,0][thefile["Flag"][range2,0]==1], bins=100, range=[0,1], alpha=0.5, color="magenta", edgecolor='magenta')
    plt.hist(ll[:,0][thefile["Flag"][range2,0]==0], bins=100, range=[0,1], alpha=0.5, color="gold", edgecolor="gold")
    plt.savefig(name_String+"plot_TrueFalse.png")
    plt.clf() 

    #roc curve to be added--> take roc from traininig
    
    fpr, tpr, _ = roc_curve(thefile["Flag"][range2,0]==1, ll[:,0])    
    print "ROC AUC:  ",auc(fpr, tpr)
    
    plt.semilogy(tpr, fpr, label='ROC curve')       
    plt.semilogy([0, 1], [0, 1], 'k--')
    plt.xticks( [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] )
    plt.grid(True, which='both')
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('broken cluster efficiency')
    plt.ylabel('broken cluster mistag rate')
    plt.title('Roc curve, auc '+str(round(auc(fpr, tpr),2))) 
    #out_file.close()
    plt.savefig(name_String+"rocCurve_TrueFalse.png")
    plt.clf()
    
    
def evaluate_onlyHIT(keras_model, npz_file, epoch):
    
    plt.clf()
    
    name_String="Eval_epoch"+str(epoch)+"_"
    
    print "\n"
    
    thefile=numpy.load(npz_file)
    sc = StandardScaler().fit(thefile["track"])
    print sc.mean_

    sc2 = StandardScaler().fit(thefile["track1"])
    print sc2.mean_
    print sc2.scale_
    
    len2=len(thefile["track"])

    range1=range(0,len2,2)
    range2=range(1,len2,2)
    
    ll=keras_model.predict([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]])
    
    
    #variabili del cluster
    ll2=thefile["hit_pos"][range2,1]
    delta=thefile["delta_TH"][range2,1]
    track_pos=thefile["track_pos"][range2,1]
    clucenter=thefile["local_cluster"][range2,1]
    cluBroken=thefile["local_brokencluster"][range2,1]
    broekn=thefile["Flag"][range2,0]
    
    print ll.shape,  ll2.shape
    
    for i in range(100):
        print ll[i], "vs", ll2[i]
        #print thefile["truth1"][i,1]
        

    plt.hist(ll[:,0], bins=100, range=[-10,10], alpha=0.5)
    plt.hist(ll2, bins=100, range=[-10,10], alpha=0.5)
    plt.savefig(name_String+"plot.png")
    plt.clf()

    plt.hist(ll[:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_EstimateOK.png")
    plt.clf() 


    plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.hist(delta, bins=100, range=[-2,2], alpha=0.5)
    plt.savefig(name_String+"plot_resolutionTrack1.png")
    plt.clf() 

    plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(delta, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_resolutionTrack.png")
    plt.clf() 

    plt.hist2d(ll2, ll[:,0],bins=100, norm=LogNorm())
    cb=plt.colorbar()
    plt.savefig(name_String+"plot_corrRainbowLog.png")
    cb.remove()
    plt.clf()


    plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(ll[:,0]-clucenter, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.hist(ll[:,0]-cluBroken, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    plt.savefig(name_String+"plot_resolutionCluster.png")
    plt.clf()


    plt.hist(ll[:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken*(broekn==1)+clucenter*(broekn==0)-ll2, bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos.png")
    plt.clf()

    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyBroken.png")
    plt.clf()

    plt.hist(ll[broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.hist(track_pos[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyFull.png")
    plt.clf()
    
    
    plt.hist(ll[:,0][broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='royalblue', label="full clusters")
    plt.hist(ll[:,0][broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='gold', label="broken clusters")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_Estimate_sovr.png")
    plt.clf() 
    
    plt.hist([ll[:,0][broekn==0]-ll2[broekn==0],ll[:,0][broekn==1]-ll2[broekn==1]], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color=['royalblue','gold'], label=["full clusters","broken clusters"], stacked=True)   
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_Estimate_stack.png")
    plt.clf()
    
    plt.hist(ll[broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlycluster.png")
    plt.clf()
    
    plt.hist(ll[broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-0.5,0.5], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-0.5,0.5], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlycluster_small.png")
    plt.clf()
    
    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyclusterBroken.png")
    plt.clf()
    
    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track position Y")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyTraBroken.png")
    plt.clf()
    

    
def Evaluate_hitOnly_wLegend(keras_model_string, npz_file, epoch):
    
    #from keras.models import Model, load_model
    
    #keras_model=load_model(keras_model_string)
    
    plt.clf()
    
    name_String="Evaluation2_"+str(epoch)+"_"
    
    print "\n"
    
    thefile=numpy.load(npz_file)
    sc = StandardScaler().fit(thefile["track"])
    print sc.mean_

    sc2 = StandardScaler().fit(thefile["track1"])
    print sc2.mean_
    print sc2.scale_
    
    len2=len(thefile["track"])

    range1=range(0,len2,2)
    range2=range(1,len2,2)
    
    #ll=keras_model.predict([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]])
    
    
    #variabili del cluster
    ll2=thefile["hit_pos"][range1,1]
    delta=thefile["delta_TH"][range1,1]
    track_pos=thefile["track_pos"][range1,1]
    clucenter=thefile["local_cluster"][range1,1]
    cluBroken=thefile["local_brokencluster"][range1,1]
    broekn=thefile["Flag"][range1,0]
    
    #print ll.shape,  ll2.shape
    
    #for i in range(100):
        #print ll[i], "vs", ll2[i]
        ##print thefile["truth1"][i,1]
        
        
    plt.hist(ll2-clucenter, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')    
    plt.savefig(name_String+"plot.png")
    plt.clf()    
    
    plt.hist(ll2[broekn==0]-clucenter[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none')    
    plt.savefig(name_String+"plot1.png")
    plt.clf()    
    
    plt.hist(ll2[broekn==1]-clucenter[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none')    
    plt.savefig(name_String+"plot2.png")
    plt.clf()   
    
    plt.hist(ll2[broekn==1]-cluBroken[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none')    
    plt.savefig(name_String+"plot3.png")
    plt.clf()   


    plt.hist([ll2[broekn==1]-cluBroken[broekn==1],ll2[broekn==0]-clucenter[broekn==0]], bins=100, range=[-2,2], alpha=0.5, edgecolor='green', stacked=True, label=['broken clusters','full clusters'])   
    plt.xlabel("hit_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot4.png")
    plt.clf()  
    
    
    plt.hist(ll2[broekn==1]-cluBroken[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="broken clusters")
    plt.hist(ll2[broekn==0]-clucenter[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="full clusters") 
    
    plt.xlabel("hit_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    
    plt.legend(loc=1)
    
    plt.savefig(name_String+"plot5.png")
    plt.clf() 
    
    plt.hist(ll2[broekn==1]-cluBroken[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='gold', label="broken clusters")
    #plt.hist(ll2[broekn==0]-clucenter[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='white', label="full clusters") 
    
    plt.xlabel("hit_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    plt.ylim(0,18000)
    
    #plt.legend(loc=1)
    
    plt.savefig(name_String+"plot5WG.png")
    plt.clf() 
    
    
    plt.hist(ll2[broekn==1]-clucenter[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='gold', label="broken clusters")
    #plt.hist(ll2[broekn==0]-clucenter[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='white', label="full clusters") 
    
    plt.xlabel("hit_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    plt.ylim(0,18000)
    
    #plt.legend(loc=1)
    
    plt.savefig(name_String+"plot5WGs.png")
    plt.clf() 
    
    
    
    plt.hist(ll2[broekn==1]-cluBroken[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='white', label="broken clusters")
    plt.hist(ll2[broekn==0]-clucenter[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', color='gold', label="full clusters") 
    
    plt.xlabel("hit_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    
    #plt.legend(loc=1)
    
    plt.savefig(name_String+"plot5WG1.png")
    plt.clf() 
    
    
    plt.hist([ll2[broekn==1]-track_pos[broekn==1],ll2[broekn==0]-track_pos[broekn==0]], bins=100, range=[-2,2], alpha=0.5, edgecolor='green', stacked=True, label=['broken clusters','full clusters'])   
    plt.xlabel("track_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot6.png")
    plt.clf()  
    
    
    plt.hist(ll2[broekn==1]-track_pos[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="broken clusters")
    plt.hist(ll2[broekn==0]-track_pos[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="full clusters") 
    
    plt.xlabel("track_position_Y-cluster_center_Y [local units]")
    plt.ylabel("# clusters")
    
    plt.legend(loc=1)
    
    plt.savefig(name_String+"plot7.png")
    plt.clf() 
    
    
    plt.hist([ll2[broekn==1]-track_pos[broekn==1],ll2[broekn==0]-track_pos[broekn==0]], bins=100, range=[-2,2], alpha=0.5, edgecolor='green', stacked=True, label=['broken clusters','full clusters'])   
    plt.xlabel("hit_position_Y-track_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot8.png")
    plt.clf()  
    
    
    plt.hist(ll2[broekn==1]-track_pos[broekn==1], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="broken clusters")
    plt.hist(ll2[broekn==0]-track_pos[broekn==0], bins=100, range=[-2,2], alpha=0.5, edgecolor='none', label="full clusters") 
    
    plt.xlabel("hit_position_Y-track_position_Y [local units]")
    plt.ylabel("# clusters")
    
    plt.legend(loc=1)
    
    plt.savefig(name_String+"plot9.png")
    plt.clf() 
    
    
    
    
    #plt.hist(ll[:,0], bins=100, range=[-10,10], alpha=0.5)
    #plt.hist(ll2, bins=100, range=[-10,10], alpha=0.5)
    #plt.savefig(name_String+"plot.png")
    #plt.clf()

    #plt.hist(ll[:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.savefig(name_String+"plot_EstimateOK.png")
    #plt.clf() 


    #plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    ##plt.hist(delta, bins=100, range=[-2,2], alpha=0.5)
    #plt.savefig(name_String+"plot_resolutionTrack1.png")
    #plt.clf() 

    #plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.hist(delta, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.savefig(name_String+"plot_resolutionTrack.png")
    #plt.clf() 

    #plt.hist2d(ll2, ll[:,0],bins=100, norm=LogNorm())
    #cb=plt.colorbar()
    #plt.savefig(name_String+"plot_corrRainbowLog.png")
    #cb.remove()
    #plt.clf()


    #plt.hist(ll[:,0]-track_pos, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.hist(ll[:,0]-clucenter, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.hist(ll[:,0]-cluBroken, bins=100, range=[-2,2], alpha=0.5, edgecolor='none')
    #plt.savefig(name_String+"plot_resolutionCluster.png")
    #plt.clf()


    #plt.hist(ll[:,0]-ll2, bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    #plt.hist(track_pos-ll2, bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    #plt.hist(clucenter-ll2, bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    #plt.hist(cluBroken*(broekn==1)+clucenter*(broekn==0)-ll2, bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    #plt.legend(loc=1)
    #plt.savefig(name_String+"plot_resolution_hitPos.png")
    #plt.clf()

    #plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    #plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    #plt.hist(clucenter[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    #plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    #plt.legend(loc=1)
    #plt.savefig(name_String+"plot_resolution_hitPos_onlyBroken.png")
    #plt.clf()

    #plt.hist(ll[broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    #plt.hist(track_pos[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track pos")
    #plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    #plt.hist(clucenter[broekn==0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='aqua',edgecolor='none', label="center w bias")
    #plt.legend(loc=1)
    #plt.savefig(name_String+"plot_resolution_hitPos_onlyFull.png")
    #plt.clf()

    
def Evaluate_hitOnly_wLegend2(keras_model_string, npz_file, epoch):
    
    from keras.models import Model, load_model
    
    keras_model=load_model(keras_model_string)
    
    plt.clf()
    
    name_String="Evaluation_bis_"+str(epoch)+"_"
    
    print "\n"
    
    thefile=numpy.load(npz_file)
    sc = StandardScaler().fit(thefile["track"])
    print sc.mean_

    sc2 = StandardScaler().fit(thefile["track1"])
    print sc2.mean_
    print sc2.scale_
    
    len2=len(thefile["track"])

    range1=range(0,len2,2)
    range2=range(1,len2,2)
    
    ll=keras_model.predict([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]])
    
    
    #variabili del cluster
    ll2=thefile["hit_pos"][range2,1]
    delta=thefile["delta_TH"][range2,1]
    track_pos=thefile["track_pos"][range2,1]
    clucenter=thefile["local_cluster"][range2,1]
    cluBroken=thefile["local_brokencluster"][range2,1]
    broekn=thefile["Flag"][range2,0]
    
    print ll.shape,  ll2.shape
    
    for i in range(100):
        print ll[i], "vs", ll2[i]
        #print thefile["truth1"][i,1]
        
    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.savefig(name_String+"plot.png")
    plt.clf() 
    
    plt.hist(ll[broekn==0][:,0]-ll2[broekn==0], bins=100, range=[-2,2], alpha=0.5, color='red', edgecolor='red', label="hit pos estimate")
    plt.savefig(name_String+"plot1.png")
    plt.clf()
    
    plt.hist([ll[broekn==1][:,0]-ll2[broekn==1],ll[broekn==0][:,0]-ll2[broekn==0]], bins=100, range=[-2,2], alpha=0.5, color=['red', 'gold'], edgecolor='red', label=['broken clusters','full clusters'], stacked=True)
    plt.savefig(name_String+"plot2.png")
    plt.clf() 

    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(track_pos[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='lime',edgecolor='none', label="track position Y")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyTraBroken.png")
    plt.clf()
    
    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(cluBroken[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='gold',edgecolor='none', label="cluster center")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyclusterBroken.png")
    plt.clf()
    
    plt.hist(ll[broekn==1][:,0]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='royalblue', edgecolor='royalblue', label="hit pos estimate")
    plt.hist(clucenter[broekn==1]-ll2[broekn==1], bins=100, range=[-2,2], alpha=0.5, color='darkorange',edgecolor='none', label="cluster center (MC)")
    plt.xlabel("estimate of hit_position_Y - hit_position_Y [local units]")
    plt.ylabel("# clusters")
    plt.legend(loc=1)
    plt.savefig(name_String+"plot_resolution_hitPos_onlyclusterBroken22.png")
    plt.clf()
          

