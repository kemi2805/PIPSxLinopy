#ifndef LINOPY
#define LINOPY

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>




class Linopy {
private:
    // Linopy stores the Matrix in this type of Vector Files. Every Matrix has 3 corresponding Files
    std::unique_ptr<std::vector<int>> A_col;
    std::unique_ptr<std::vector<int>> A_row;
    std::unique_ptr<std::vector<double>> A_data;

    std::unique_ptr<std::vector<int>> B_col;
    std::unique_ptr<std::vector<int>> B_row;
    std::unique_ptr<std::vector<double>> B_data;

    std::unique_ptr<std::vector<int>> BL_col;
    std::unique_ptr<std::vector<int>> BL_row;
    std::unique_ptr<std::vector<double>> BL_data;

    std::unique_ptr<std::vector<int>> C_col;
    std::unique_ptr<std::vector<int>> C_row;
    std::unique_ptr<std::vector<double>> C_data;

    std::unique_ptr<std::vector<int>> D_col;
    std::unique_ptr<std::vector<int>> D_row;
    std::unique_ptr<std::vector<double>> D_data;

    std::unique_ptr<std::vector<int>> DL_col;
    std::unique_ptr<std::vector<int>> DL_row;
    std::unique_ptr<std::vector<double>> DL_data;

    //Linopy stores the constraints and the objective function in vectors. 
    std::unique_ptr<std::vector<double>> c;  // Objective function constraints
    std::unique_ptr<std::vector<double>> b;  // Equality Constraints
    std::unique_ptr<std::vector<double>> bl; // Equality Linking Constraints

    std::unique_ptr<std::vector<double>> dl; // All kind of inequality Constraints
    std::unique_ptr<std::vector<double>> idl; //Indikator Vektor für dl. Nicht von linopy, für PIPS
    std::unique_ptr<std::vector<double>> du;
    std::unique_ptr<std::vector<double>> idu; //Indikator Vektor für du. Nicht von linopy, für PIPS
    std::unique_ptr<std::vector<double>> dll;
    std::unique_ptr<std::vector<double>> idll; //Indikator ...
    std::unique_ptr<std::vector<double>> dlu;
    std::unique_ptr<std::vector<double>> idlu; //Indikator ..

    std::unique_ptr<std::vector<double>> xl;
    std::unique_ptr<std::vector<double>> ixl; //Indikator ...
    std::unique_ptr<std::vector<double>> xu;
    std::unique_ptr<std::vector<double>> ixu; //Indikator ...

    std::unique_ptr<std::vector<int>> xvec; // The Column numbers each block uses. A is completly in Block 0
    std::unique_ptr<std::vector<int>> xvec_A; // Das xvec aber nur für A


    int nScenario;
    int id = 0; // Tracks the id of each block. First get 0, the children get from 1 to nScenarios
    std::string PathToFile; // Exact Name and Block of the File
    std::string Filepath; // General Directory of the Files

    void Set_equality_Vector(); // Those take the Data and Set every variable
    void Set_inequality_Vector();
    void Set_Matrix();


    void SetPathToFile(int id, std::string name);

    
    std::pair<std::unique_ptr<std::vector<double>>,std::unique_ptr<std::vector<double>>> GetDataFromFileInequalityConstraints(std::string name, bool linking = false);

    template<typename T>
    std::unique_ptr<std::vector<T>> GetDataFromFileEqualityConstraints(std::string name, bool linking = false);
    
    std::unique_ptr<std::vector<int>> GetDataFromRowFile(std::string name, bool linking = false);
    std::unique_ptr<std::vector<int>> GetDataFromColFile(std::string name, bool linking = false);
    std::unique_ptr<std::vector<double>> GetDataFromDataFile(std::string name, bool linking = false);
    void GetDataForxvec(std::string name = "x", bool linking = false);


public:
    Linopy(std::string filepath, int local_id, int id_size = 0);
    Linopy();
    //~Linopy(); Bekomme abgrundtief viele Fehler.

    int GetId();
    void SetId(int Id_number);
    int Get_n() const;
    int Get_my() const;
    int Get_mz() const;
    int Get_myl() const;
    int Get_mzl() const;

    int Get_nnz_A() const;
    int Get_nnz_B() const;
    int Get_nnz_C() const;
    int Get_nnz_D() const;
    int Get_nnz_BL() const;
    int Get_nnz_DL() const;

    //Linopy_Init greift auf die Funktionen zu. Auf die oben auch, aber die sind hauptsächlich dafür geschrieben.
    // They get the Pointer. The will then be moved to PIPS
    std::unique_ptr<std::vector<double>> Get_equality_vector() ;
    std::unique_ptr<std::vector<double>> Get_equality_linking_vector() ;
    std::unique_ptr<std::vector<double>> Get_objective_function() ;

    std::unique_ptr<std::vector<double>> Get_upper_inequality_constraint_vector() ; // du
    std::unique_ptr<std::vector<double>> Get_lower_inequality_constraint_vector() ;
    std::unique_ptr<std::vector<double>> Get_upper_inequality_constraint_linking_vector() ;
    std::unique_ptr<std::vector<double>> Get_lower_inequality_constraint_linking_vector() ;

    std::unique_ptr<std::vector<double>> Get_upper_inequality_constraint_vector_indicator() ; // du
    std::unique_ptr<std::vector<double>> Get_lower_inequality_constraint_vector_indicator() ;
    std::unique_ptr<std::vector<double>> Get_upper_inequality_constraint_linking_vector_indicator() ;
    std::unique_ptr<std::vector<double>> Get_lower_inequality_constraint_linking_vector_indicator() ;

    std::unique_ptr<std::vector<double>> Get_upper_inequality_vector() ;
    std::unique_ptr<std::vector<double>> Get_lower_inequality_vector() ;

    std::unique_ptr<std::vector<double>> Get_upper_inequality_vector_indicator() ;
    std::unique_ptr<std::vector<double>> Get_lower_inequality_vector_indicator() ;

    std::unique_ptr<std::vector<int>> Get_A_row() ;
    std::unique_ptr<std::vector<int>> Get_A_col() ;
    std::unique_ptr<std::vector<double>> Get_A_data() ;

    std::unique_ptr<std::vector<int>> Get_B_row() ;
    std::unique_ptr<std::vector<int>> Get_B_col() ;
    std::unique_ptr<std::vector<double>> Get_B_data() ;

    std::unique_ptr<std::vector<int>> Get_BL_row() ;
    std::unique_ptr<std::vector<int>> Get_BL_col() ;
    std::unique_ptr<std::vector<double>> Get_BL_data() ;

    std::unique_ptr<std::vector<int>> Get_C_row();
    std::unique_ptr<std::vector<int>> Get_C_col() ;
    std::unique_ptr<std::vector<double>> Get_C_data() ;

    std::unique_ptr<std::vector<int>> Get_D_row() ;
    std::unique_ptr<std::vector<int>> Get_D_col() ;
    std::unique_ptr<std::vector<double>> Get_D_data() ;

    std::unique_ptr<std::vector<int>> Get_DL_row() ;
    std::unique_ptr<std::vector<int>> Get_DL_col() ;
    std::unique_ptr<std::vector<double>> Get_DL_data() ;

    int* Get_xvec_A();
    int Get_xvec_A_size();
    void Set_xvec_A(int* xvec_A0, int size);
    std::unique_ptr<std::vector<int>> Get_xvec();

    std::unique_ptr<std::vector<int>> Combine_yvec(std::unique_ptr<std::vector<int>>yvec1,std::unique_ptr<std::vector<int>>yvec2);

    std::unique_ptr<std::vector<int>> Transform_Matrix_Rows(std::unique_ptr<std::vector<int>> M_Row);//Those do what they say
    void Transform_Matrix_Cols();
    void Transform_Matrix_Cols_A(std::vector<int>::iterator Col_begin, std::vector<int>::iterator Col_end);
    void Transform_Matrix_Cols_B(std::vector<int>::iterator Col_begin, std::vector<int>::iterator Col_end);








};

#endif