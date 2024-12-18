############
PhosphoLingo
############

Introduction
############
PhosphoLingo is an advanced phosphorylation predictor using deep learning and Protein Language Models. It predicts likely phosphorylation sites from sequence alone. Given the used data format, it can easily be extended to other post-translational modifications (PTMs).

Supported functionality
***********************
The PhosphoLingo tool can be called for three different goals:

- training a new model (**train**) and evaluating it on a test set (included)
- make predictions for new sequences using an existing model (**predict**)
- visualize important features using an existing model (**visualize**)

Input data format
#################
The input sequences should be provided using a FASTA format, with the positive and negative phosphorylation (or by extension, any PTM) sites indicated. This can be done by adding a trailing ``#`` to the modified residue for positive annotations, or a ``@`` for negative annotations.

Dummy example for Y phosphorylation, with two positive annotations, four negative annotations, and three Y residues without annotations:

.. code-block::

    >dummy_protein
    MKPWETDY#MGPFRKIY@Y@WIVFYESGSDMDCNY#CYNFKGEFITQSICHPNDPSKEGDQARTWCIS
    ETLRSRTTFYDLSIGLTKSFFFNRY@CIWFCFIWLSDDRMTKY@

When making predictions/visualizing for protein sequences of which no annotations might be available, the sites to be predicted should be indicated by one of the two tokens (``#`` or ``@``) - in this scenario, labels are ignored so both will work equally.

Configuration file format
#########################
.. _configs: https://github.com/jasperzuallaert/PhosphoLingo/tree/master/configs
The architecture and training hyperparameters, dataset locations, and more information is stored in a configuration file, in a ``.json`` format. Example configurations can be found in the configs_ directory. The different options are given below.

====================== ========================== ===
name                   default value
====================== ========================== ===
training_set           "Scop3P/ST/PF"             The location of FASTA files to train, validate and test (default) or to train and validate (if ``test_set != "default"``)
test_set               "default"                  "default" if regular training/validation/test scheme is followed, otherwise the location of the ``fold[0-9].fasta`` files for evaluation
test_fold              0                          Only used if ``test_set != "default"``. Chooses the ``fold[0-9].fasta`` file to use for evaluation.
save_model             False                      If True, the model will be kept after training.
representation         "ProtTransT5_XL_UniRef50"  The representation to use. Supported: ``onehot``, ``ESM1_small``, ``ESM1b``, ``ESM2_150M``, ``ESM2_650M``, ``ESM2_3B``, ``CARP_640M``, ``ProtTransT5_XL_UniRef50``, ``Ankh_base``, ``Ankh_large``
freeze_representation  True                       Whether to freeze the representation during training or fine-tune it. All models in the paper are trained with this parameter set to ``True``
receptive_field        65                         The number of residues (centered around the candidate P-site) that are given to the CNN
conv_depth             3                          The number of convolutional blocks in the architecture
conv_width             9                          The filter size of the convolutional layers
conv_channels          200                        The number of output channels (= number of filters) per convolutional layer
final_fc_neurons       64                         The number of neurons in the final fully connected layer
dropout                0.2                        Only used if ``batch_norm == False``. Dropout *p* for all dropout layers
max_pool_size          1                          The pooling size used in the pooling layers
batch_norm             False                      Chooses whether to use batch normalization as regularization (if ``True``), or dropout (if ``False``)
batch_size             4                          The number of protein fragments in each batch
warm_up_epochs         1.5                        The number of epochs for learning rate warm-up (linear from 0.1*LR)
learning_rate          1e-4                       Learning rate for AdamW

pos_weight             1                          The weight for positive samples in the loss calculation
max_epochs             10                         The maximum number of epochs, if no early stopping occurs
====================== ========================== ===

Usage
*****
For hints on how to run PhosphoLingo:

.. code-block::

    python phospholingo -h

.. code-block::

    usage: phospholingo [-h] {train,predict,visualize} ...

    positional arguments:
      {train,predict,visualize}
                            the desired PhosphoLingo subprogram to run
        train               train a new prediction model
        predict             predict using an existing model
        visualize           calculate SHAP values using an existing model

    optional arguments:
      -h, --help            show this help message and exit

Train a new model
#################
Usage:

.. code-block ::

    usage: phospholingo train [-h] json

    positional arguments:
      json        the .json configuration file for training a new model

    optional arguments:
      -h, --help  show this help message and exit

To train a new model, supply a ``.json`` file with the desired configuration. The AUPRC, AUROC, and precisions at recall of 0.8 and 0.6 will be logged in the resulting directory. If specified in the configuration file, the checkpoint of the model will also be saved.

Training can be done via two data setups:

- **(default)** training/validation/test sets: The default training run. This is achieved by setting ``test_set`` to ``default`` in the config. In this case, training will be done on the ``train.fasta`` file in the specified data directory (``dataset`` in the config), early stopping will be done using ``valid.fasta``, and test metrics are computed on the ``test.fasta`` data.
- cross-dataset evaluation: Specifically to reproduce results in the paper or to check model transferability between datasets. This is achieved by setting ``test_set`` to the desired data directory on which evaluation should be done. Additionally, specify the fold to test on by setting ``test_fold`` to any number between 0 and 9. The ``fold[0-9].fasta`` file will be used for evaluation, and all proteins present will be removed from the training and validation sets.



Predict using an existing model
###############################
Usage:

.. code-block ::

    usage: phospholingo predict [-h] model dataset out

    positional arguments:
      model       the location of the saved model
      dataset     the dataset for which to make predictions
      out         the output file, will be written in a csv format

    optional arguments:
      -h, --help  show this help message and exit

You can make predictions on an unseen ``dataset``, using a pretrained prediction ``model``, and writing results to an ``out`` csv file. As indicated before, the dataset should be in a FASTA format, and sites to be predicted should be followed by either a ``#`` or ``@`` symbol. The actual annotations are ignored, so either symbol will work equivalently.


Mention data format again

Visualize important features using an existing model
####################################################
Usage:

.. code-block ::

    usage: phospholingo visualize [-h] model dataset out_values out_img

    positional arguments:
      model       the location of the saved model
      dataset     the dataset for which to visualize important features
      out_values  the output SHAP scores file, will be written in a txt format
      out_img     the normalized average SHAP scores per position, as an image file

    optional arguments:
      -h, --help  show this help message and exit

You can make visualizations for a pretrained prediction ``model``, on an unseen ``dataset``. Output values will be stored in ``out_values`` (.txt format), and an image will be generated to ``out_img`` (.jpg, .png, .svg, ...).

Setting the maximum system batch size
*************************************
.. _utils: https://github.com/jasperzuallaert/PhosphoLingo/blob/master/phospholingo/utils.py
As Protein Language Models can be very resource-heavy to use, especially when considering the larger models and when also fine-tuning them during training, the user can set their maximum batch size for specific situations. This is done in the ``get_gpu_max_batchsize`` function in utils_. Users can redefine this function so that appropriate batch sizes are returned for their system. A non-optimized example for different batch sizes using different representations is implemented, though this has not been thoroughly optimized.


Extra files
***********
Pre-trained phosphorylation models (``.ckpt`` format) can be downloaded from following locations. The models are trained on the Scop3P-ST-PF and Scop3P-Y-PF datasets.

====================== ======= =========
Model                  Targets Link
====================== ======= =========
ProtT5-XL-U50           ST     https://huggingface.co/JasperZ/PhosphoLingo_ST/resolve/main/PhosphoLingo_ST_new.ckpt
ProtT5-XL-U50           Y      https://huggingface.co/JasperZ/PhosphoLingo_Y/resolve/main/PhosphoLingo_Y_new.ckpt
====================== ======= =========

.. _data: https://github.com/jasperzuallaert/PhosphoLingo/blob/master/data/
.. _predictions: https://github.com/jasperzuallaert/PhosphoLingo/blob/master/predictions/
Datasets (FASTA format with ``#`` and ``@`` annotations) used in this study are found in data_.

Configuration files (``.json`` format) can be found in configs_. If you want to run the preset configurations, you should only change the following parameters: ``training_set``, ``test_set``, ``test_fold``, and ``save_model``.

Prediction files for the human proteome as found in SwissProt (version 11/2022) are found in predictions_.

Cite
****
UNDER CONSTRUCTION
