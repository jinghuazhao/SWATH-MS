# 1. Load appropriate modules and environment
ssh login-gpu.hpc.cam.ac.uk
module purge
module load rhel7/default-gpu
module unload cuda/8.0
module load python/3.6 cuda/10.0 cudnn/7.5_cuda-10.0

# 2. Install latest Tensorflow in a python virtualenv
python -m venv $HOME/.conda/tensorflow-gpu
source $HOME/.conda/tensorflow-gpu/bin/activate
pip install tensorflow-gpu

# 3. Test your installation (be in a interative session
# via sintr or submit a job to a gpu compute node)
cat << 'EOF' > helloworld.py
#!/usr/bin/env python
import tensorflow as tf
mnist = tf.keras.datasets.mnist

(x_train, y_train),(x_test, y_test) = mnist.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0

model = tf.keras.models.Sequential([
          tf.keras.layers.Flatten(input_shape=(28, 28)),
          tf.keras.layers.Dense(512, activation=tf.nn.relu),
          tf.keras.layers.Dropout(0.2),
          tf.keras.layers.Dense(10, activation=tf.nn.softmax)
          ])
model.compile(optimizer='adam',
          loss='sparse_categorical_crossentropy',
          metrics=['accuracy'])

model.fit(x_train, y_train, epochs=5)
model.evaluate(x_test, y_test)
EOF
