#include "readCSV.hpp"


std::vector<std::vector<std::string>> readCsvFile(std::string filename)
{
    std::ifstream       file;
    file.open(filename, std::ifstream::in);
    if (file.is_open()) {
        std::cout<<"file Open" <<std::endl;

    }
    else {
        // show message:
        std::cout << "Error opening file";
    }


    std::vector<std::vector<std::string>> read;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
        std::vector<std::string> vec;//vector of 8 strings

        for ( int i= 0; i < (*loop).size() ; i++)
        {
            vec.push_back((*loop)[i]);
           // std::cout << (*loop)[i] ;
        }
        read.push_back(vec);

    }
    return read;
}


std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}
