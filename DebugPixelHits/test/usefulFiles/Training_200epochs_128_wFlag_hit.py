import numpy
#import numpy
seed = 7 
numpy.random.seed(seed)

import sys

from sklearn.preprocessing import StandardScaler
from evaluation_function import evaluate, evaluateONLY_TrueFalse, evaluate_onlyHIT

name_String=sys.argv[1].split('.')[0]+"Eval_wFlag_"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

thefile=numpy.load(sys.argv[1])
sc = StandardScaler().fit(thefile["track"])
print sc.mean_

sc2 = StandardScaler().fit(thefile["track1"])
print sc2.mean_
print sc2.scale_

len2=len(thefile["track"])
len1=len(thefile["track"])/2


range1=range(0,len2,2)
range2=range(1,len2,2)

#plt.hist(thefile["hit_pos"][range1,1], bins=100, range=[-10,10])
#plt.savefig(name_String+"truth_train.png")

#plt.clf()
#plt.hist(thefile["hit_pos"][range2,1], bins=100, range=[-10,10])
#plt.savefig(name_String+"truth_val.png")

plt.clf()

from keras.layers import Dense, Flatten, Convolution2D, MaxPooling2D, merge
from keras.models import Model, load_model
from keras.layers import Input 
from random import randint 
import keras

class wHistory(keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs={}): 
        if (epoch%5==4):
            self.model.save(sys.argv[1].split('.')[0]+"weights_128_wFlag_at_epoch_hit"+str(epoch)+".h5") 
            evaluate_onlyHIT(self.model, sys.argv[1], epoch)
        #string_name="eval3Test_at_epoch"+str(epoch)+".h5"
        
wHistory_ = wHistory()
reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=5, min_lr=0.001)

Inputs=[Input(shape=(3,)) ]
Inputs+=[Input(shape=(7, ))]
Inputs+=[Input(shape=(21, 7, 2)) ]

image=Inputs[2]

image  = Convolution2D(4, 3, 3, init='lecun_uniform',  activation='relu')(image)
image  = Convolution2D(4, 2, 2, init='lecun_uniform',  activation='relu')(image)
#image  = Convolution2D(64, 3, 3, init='lecun_uniform',  activation='relu')(image)
image  = MaxPooling2D(pool_size=(2, 2))(image)
image = Flatten()(image)

x = merge( [Inputs[0] , Inputs[1], image] , mode='concat')


x=  Dense(100, activation='relu',init='lecun_uniform')(x)
x=  Dense(50, activation='relu',init='lecun_uniform')(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
x=  Dense(20, activation='relu',init='lecun_uniform')(x)

predictions = Dense(1,init='lecun_uniform', name='out1')(x)
#predictions2 = Dense(1,activation='sigmoid',init='lecun_uniform', name='out2')(x)
model = Model(input=Inputs, output=predictions)

model.compile(loss='mse', optimizer='adam',metrics=['mae', 'mse'])

model.summary()

history=model.fit([sc.transform(thefile["track"][range1]), sc2.transform(thefile["track1"][range1]), thefile["image"][range1]], (thefile["hit_pos"][range1,1]) ,  nb_epoch=200, verbose=1,batch_size=128, 
                  callbacks=[wHistory_, reduce_lr], validation_data=([sc.transform(thefile["track"][range2]), sc2.transform(thefile["track1"][range2]), thefile["image"][range2]], (thefile["hit_pos"][range2,1])))

# list all data in history
print(history.history.keys())
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.savefig(name_String+"plot_history.png")
plt.clf()


plt.semilogy(history.history['loss'])
plt.semilogy(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.savefig(name_String+"plot_history_log.png")
plt.clf()


