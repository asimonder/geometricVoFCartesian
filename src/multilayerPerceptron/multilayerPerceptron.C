/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------

License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multilayerPerceptron.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::multilayerPerceptron::multilayerPerceptron(std::string filename):
  fileName_(filename)
{
  std::string fName=fileName_+"_activations.txt";
  std::ifstream file(fName);
  if (!file.is_open())
    std::cout << "Error opening file "<<fName<<" \n";
  std::string line;
  while (getline(file, line))
    {
      std::stringstream ss (line);
      activations_.push_back(ss.str());
    }
  file.close();
  
  NLayers_=activations_.size();
  
  for(int i=0;i<NLayers_;i++)
    {
      std::ostringstream oss;
      oss<<fileName_<<"_weights"<<i<<".txt";
      std::ifstream fileW(oss.str());
      if (!fileW.is_open())
	std::cout << "Error opening file "<<oss.str()<<" \n";
      std::vector<std::vector<double>> layer;
      while (getline(fileW, line))
	{
	  double itmp;
	  std::vector<double> vec;
	  std::stringstream ss (line);
	  while (ss >> itmp)
	    {
	      std::string tmpstr;
	      vec.push_back(itmp);
	      if (!getline (ss, tmpstr,','))     //read to ',' w/tmpstr 
		break;                          //done if no more ',' 
	    }
	  layer.push_back(vec);
	}
      weights_.push_back(layer);
      fileW.close();

      std::ostringstream oss2;
      oss2<<fileName_<<"_biases"<<i<<".txt";
      std::ifstream fileB(oss2.str());
      if (!fileB.is_open())
	std::cout << "Error opening file "<<oss2.str()<<" \n";
      std::vector<double> vec;
      while (getline(fileB, line))
	{
	  double itmp;
	  std::stringstream ss (line);
	  ss >> itmp;
	  vec.push_back(itmp);
	}
      biases_.push_back(vec);
      fileB.close();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
std::vector<std::string> Foam::multilayerPerceptron::activations()
{
  return activations_;
}

std::vector<std::vector<std::vector<double>>> Foam::multilayerPerceptron::weights()
{
  return weights_;
}

std::vector<std::vector<double>> Foam::multilayerPerceptron::biases()
{
  return biases_;
}

std::vector<double> Foam::multilayerPerceptron::predict(const std::vector<double>& data)
{
  //std::cout<<"predicting...\n";
  std::vector<double> inputLayer(data);
  //inputLayer.push_back(data);
  for(int i=0;i<NLayers_;i++)
    {
      std::vector<double> outputLayer;
      int NNeurons=weights_[i].size();
      for(int j=0;j<NNeurons;j++)
	{
	  double sum=0;
	  /*for(int k=0;k<inputLayer.size();k++)
	    {
	      sum+=inputLayer[k]*weights_[i][j][k];
	      }*/
	  sum=std::inner_product(inputLayer.begin(),inputLayer.end(),weights_[i][j].begin(),0.0);
	  sum+=biases_[i][j];
	  if (activations_[i]=="tanh")
	    {
	    sum=tanh(sum);
	    }
	  else if (activations_[i]=="relu")
	    {
	    sum=max(0.0,sum);
	    }
	  else if (activations_[i]=="sigmoid")
	    {
	      //Info<<"Calculating sigmoid..."<<endl;
	      sum=1./(1.+exp(-sum));
	    }
	  else if (activations_[i]=="symsigmoid")
	    {
	      //Info<<"Calculating sigmoid..."<<endl;
	      //a = 2 ./ (1 + exp(-2*n)) - 1;
	      sum=2./(1 + exp(-2*sum)) - 1.0;
	    }

	  outputLayer.push_back(sum);
	}
      inputLayer.clear();
      inputLayer=outputLayer;
    }
  return inputLayer;
}
// ************************************************************************* //
