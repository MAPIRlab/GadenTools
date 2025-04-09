# GadenTools
Python API to access the results of GADEN simulations directly.

GadenTools is a Python module that lets you access the results of [GADEN](https://github.com/MAPIRlab/gaden) simulations directly from your Python code. It eliminates the need to install and use ROS, which up until now has been the only way to interact with the results from this simulator. The intention is for this API to create a more accessible way for researchers to leverage existing resources of high fidelity gas dispersion simulations, with a special enphasis on mobile robotics.

If you use GadenTools in your research, you can use the following BibTeX for the citation:
```
@INPROCEEDINGS{ojeda_robot2022_gadentools,
     author = {Ojeda, Pepe and Ruiz-Sarmiento, J. R. and Monroy, Javier and Gonzalez-Jimenez, Javier},
      title = {GadenTools: A Toolkit for Testing and Simulating Robotic Olfaction Tasks With Jupyter Notebook Support},
  booktitle = {Fifth Iberian Robotics Conference (ROBOT2022)},
       year = {2022},
   location = {Zaragoza, Spain},
        url = {https://mapir.uma.es/papersrepo/2022/2022_ojeda_ROBOT2022_GadenTools.pdf},
        doi = {10.1007/978-3-031-21062-4_14}
}
```

# Jupyter Notebook
GadenTools is shipped as an easy-to-use Python module, and accompanied by an extensive tutorial in the form of a [Jupyter Notebook](https://colab.research.google.com/drive/1Xj7rrsmeDa_dS3Ru_UIhhzlaifGH6GS4?usp=sharing#scrollTo=5N0zIiKtHkMY), designed to illustrate both the basic interaction with the API and some more advanced examples of how to exploit its capabilities to obtain powerful visualizations and analysis of the data.

![image](https://user-images.githubusercontent.com/5920310/177741853-5a0cadbf-f938-4ccc-bb9b-91a838386c07.png)


# Dataset
To illustrate different use cases of GadenTools, we employ a publicly availabe [Dataset](https://mapir.isa.uma.es/mapirwebsite/?p=1708) of simulated airflows and gas dispersion generated with GADEN, in detailed 3D models of real houses (made with robotics in mind).
![image](https://mapir.isa.uma.es/mapirwebsite/wp-content/uploads/2022/04/vectorglyphs.png)
