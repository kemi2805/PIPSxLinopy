#include "../../Core/Interface/PIPSIPMppInterface.hpp"
#include "../../Core/Options/PIPSIPMppOptions.h"
#include "DistributedInputTree.h"


#include "mpi.h"

//added Packages for reading directly from the files
#include <fstream>
#include <string>
#include "LinopyModel.hpp"

extern "C" typedef int (* FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (* FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (* FVEC)(void* user_data, int id, double* vec, int len);


/** Problem parameters and data */
class ProbData {
public: //data
    int nScenarios;

public: //methods
    explicit ProbData(int nScenarios) {
        this->nScenarios = nScenarios;
    };

//  ~ProbData();
};
    


    std::vector<std::unique_ptr<std::vector<long long>>> A_col;
    std::vector<std::unique_ptr<std::vector<long long>>> A_row;
    std::vector<std::unique_ptr<std::vector<double>>> A_data;

    std::vector<std::unique_ptr<std::vector<long long>>> B_col;
    std::vector<std::unique_ptr<std::vector<long long>>> B_row;
    std::vector<std::unique_ptr<std::vector<double>>> B_data;

    std::vector<std::unique_ptr<std::vector<long long>>> BL_col;
    std::vector<std::unique_ptr<std::vector<long long>>> BL_row;
    std::vector<std::unique_ptr<std::vector<double>>> BL_data;

    std::vector<std::unique_ptr<std::vector<long long>>> C_col;
    std::vector<std::unique_ptr<std::vector<long long>>> C_row;
    std::vector<std::unique_ptr<std::vector<double>>> C_data;

    std::vector<std::unique_ptr<std::vector<long long>>> D_col;
    std::vector<std::unique_ptr<std::vector<long long>>> D_row;
    std::vector<std::unique_ptr<std::vector<double>>> D_data;

    std::vector<std::unique_ptr<std::vector<long long>>> DL_col;
    std::vector<std::unique_ptr<std::vector<long long>>> DL_row;
    std::vector<std::unique_ptr<std::vector<double>>> DL_data;

    //Linopy stores the constraints and the objective function in vectors. 
    std::vector<std::unique_ptr<std::vector<double>>> c;  // Objective function constraints
    std::vector<std::unique_ptr<std::vector<double>>> b;  // Equality Constraints
    std::vector<std::unique_ptr<std::vector<double>>> bl; // Equality Linking Constraints

    std::vector<std::unique_ptr<std::vector<double>>> dl; // All kind of inequality Constraints
    std::vector<std::unique_ptr<std::vector<double>>> idl; //Indikator Vektor für dl. Nicht von linopy, für PIPS
    std::vector<std::unique_ptr<std::vector<double>>> du;
    std::vector<std::unique_ptr<std::vector<double>>> idu; //Indikator Vektor für du. Nicht von linopy, für PIPS
    std::vector<std::unique_ptr<std::vector<double>>> dll;
    std::vector<std::unique_ptr<std::vector<double>>> idll; //Indikator ...
    std::vector<std::unique_ptr<std::vector<double>>> dlu;
    std::vector<std::unique_ptr<std::vector<double>>> idlu; //Indikator ..

    std::vector<std::unique_ptr<std::vector<double>>> xl;
    std::vector<std::unique_ptr<std::vector<double>>> ixl; //Indikator ...
    std::vector<std::unique_ptr<std::vector<double>>> xu;
    std::vector<std::unique_ptr<std::vector<double>>> ixu; //Indikator ...


    std::vector<long long> n;
    std::vector<long long> my;
    std::vector<long long> myl;
    std::vector<long long> mz;
    std::vector<long long> mzl;

    std::vector<long long> nnzA;
    std::vector<long long> nnzB;
    std::vector<long long> nnzC;
    std::vector<long long> nnzD;
    std::vector<long long> nnzBL;
    std::vector<long long> nnzDL;
    int nSize(void*, int id, int* nnz) {
        *nnz = n.at(id);
        return 0;
    }
    int mySize(void*, int id, int* nnz) {
        *nnz = my.at(id);
        return 0;
    }
    int mzSize(void*, int id, int* nnz) {
        *nnz = mz.at(id);
        return 0;
    }
    int mylSize(void*, int id, int* nnz) {
        *nnz = myl.at(0);
        return 0;
    }
    int mzlSize(void*, int id, int* nnz) {
        *nnz = mzl.at(0);
        return 0;
    }
    int nnzMatEqStage1(void*, int id, int* nnz) {
        *nnz = nnzA.at(id);
        return 0;
    }
    int nnzMatIneqStage1(void*, int id, int* nnz) {
        *nnz = nnzC.at(id);
        return 0;
    }
    int nnzMatEqStage2(void*, int id, int* nnz) {
        *nnz = nnzB.at(id);
        return 0;
    }
    int nnzMatIneqStage2(void*, int id, int* nnz) {
        *nnz = nnzD.at(id);
        return 0;
    }
    int nnzMatEqLink(void*, int id, int* nnz) {
        *nnz = nnzBL.at(id);
        return 0;
    }
    int nnzMatIneqLink(void*, int id, int* nnz) {
        *nnz = nnzDL.at(id);
        return 0;
    }
    int nnzAllZero(void*, int, int* nnz) {
        *nnz = 0;
        return 0;
    }
    int vecAllZero(void*, int, double* vec, int len) {
        return 0;
    }
    int vecEqRhs(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < b.at(id)->size(); i++ ) 
            vec[i] = b.at(id)->at(i);
        return 0;
    }
    int vecEqRhsLink(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < bl.at(id)->size(); i++ ) 
            vec[i] = bl.at(0)->at(i);
        return 0;
    }
    int vecIneqRhs(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < du.at(id)->size(); i++ ) 
            vec[i] = du.at(id)->at(i);        
        return 0;
    }
    int vecIneqRhsActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < idu.at(id)->size(); i++ ) 
            vec[i] = idu.at(id)->at(i);           
        return 0;
    }
    int vecIneqLhs(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < dl.at(id)->size(); i++ ) 
            vec[i] = dl.at(id)->at(i);   
        return 0;
    }
    int vecIneqLhsActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < idl.at(id)->size(); i++ ) 
            vec[i] = idl.at(id)->at(i);   
        return 0;
    }
    int vecIneqRhsLink(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < dlu.at(id)->size(); i++ ) 
            vec[i] = dlu.at(0)->at(i);   
        return 0;
    }
    int vecIneqRhsLinkActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < idlu.at(id)->size(); i++ ) 
            vec[i] = idlu.at(0)->at(i);           
        return 0;
    }
    int vecIneqLhsLink(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < dll.at(id)->size(); i++ ) 
            vec[i] = dll.at(0)->at(i);
        return 0;
    }
    int vecIneqLhsLinkActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < idll.at(id)->size(); i++ ) 
            vec[i] = idll.at(0)->at(i);   
        return 0;
    }
    int vecObj(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < c.at(id)->size(); i++ ) 
            vec[i] = c.at(id)->at(i);   
        return 0;
    }
    int vecXlb(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < xl.at(id)->size(); i++ ) 
            vec[i] = xl.at(id)->at(i);   
        return 0;
    }
    int vecXub(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < xu.at(id)->size(); i++ ) 
            vec[i] = xu.at(id)->at(i);   
        return 0;
    }
    int vecXlbActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < ixl.at(id)->size(); i++ ) 
            vec[i] = ixl.at(id)->at(i);   
        return 0;
    }
    int vecXubActive(void*, int id, double* vec, int len) {
        for(long unsigned int i = 0; i < ixu.at(id)->size(); i++ ) 
            vec[i] = ixu.at(id)->at(i);   
        return 0;
    }
    int matAllZero(void*, int, int* krowM, int* , double*) {
        krowM[0] = 0;
        krowM[1] = 0;
        return 0;
    }
    int matEqStage1(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < A_row.at(id)->size(); i++) {
            krowM[i] = A_row.at(id)->at(i);
        }
        for(long unsigned int i = 0; i < A_col.at(id)->size(); i++) {
            jcolM[i] = A_col.at(id)->at(i);
            M[i] = A_data.at(id)->at(i);
        }
        return 0;
    }
    int matIneqLink(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < DL_row.at(id)->size(); i++) 
            krowM[i] = DL_row.at(id)->at(i);
        for(long unsigned int i = 0; i < DL_col.at(id)->size(); i++) {
            jcolM[i] = DL_col.at(id)->at(i);
            M[i] = DL_data.at(id)->at(i);
        }
        return 0;
    }
    int matIneqStage1(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < C_row.at(id)->size(); i++)
            krowM[i] = C_row.at(id)->at(i);
        for(long unsigned int i = 0; i < C_col.at(id)->size(); i++) {
            jcolM[i] = C_col.at(id)->at(i);
            M[i] = C_data.at(id)->at(i);
        }
        return 0;
    }
    int matEqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < B_row.at(id)->size(); i++)
            krowM[i] = B_row.at(id)->at(i);
        for(long unsigned int i = 0; i < B_col.at(id)->size(); i++) {
            jcolM[i] = B_col.at(id)->at(i);
            M[i] = B_data.at(id)->at(i);
        }
        return 0;
    }
    int matIneqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < D_row.at(id)->size(); i++) 
            krowM[i] = D_row.at(id)->at(i);
        for(long unsigned int i = 0; i < D_col.at(id)->size(); i++) {
            jcolM[i] = D_col.at(id)->at(i);
            M[i] = D_data.at(id)->at(i);
        }
        return 0;
    }
    int matEqLink(void*, int id, int* krowM, int* jcolM, double* M) {
        for(long unsigned int i = 0; i < BL_row.at(id)->size(); i++) 
            krowM[i] = BL_row.at(id)->at(i);
        for(long unsigned int i = 0; i < BL_col.at(id)->size(); i++) {
            jcolM[i] = BL_col.at(id)->at(i);
            M[i] = BL_data.at(id)->at(i);
        }
        return 0;
    }


int main(int argc, char** argv) {
   
    // This is for the Parallel Processor
    
   
    // This has to be the same number as the number of Blocks you put in Linopy.
    // Linopy gives an extra block out. Substract 1 from it.
    std::string nScenariostring = argv[2];
    int nScenarios = stoi(nScenariostring); //Argument aus der Main-Funktion
    std::string Filepath = argv[3];
    std::cout << "nScenario = " << nScenarios << std::endl;
    std::cout << "Filepath = " << Filepath << std::endl;

    // Hier werden alle Blöcke abgelesen und in Linopy_Storage mit der Klasse Linopy gespeichert und bearbeitet
    std::vector<Linopy> Linopy_Storage;
    for(int i = 0; i <= nScenarios; i++) {
        Linopy_Storage.push_back(Linopy(Filepath, i, nScenarios));
        std::cout << "Linopy_Storage " << i << " wurde kreirt" << std::endl;
    }

    for(int i = 0; i <= nScenarios; i++) {
        Linopy_Storage.at(i).Set_xvec_A(Linopy_Storage.at(0).Get_xvec_A(),Linopy_Storage.at(0).Get_xvec_A_size());
        Linopy_Storage.at(i).Transform_Matrix_Cols();
    }

    // Da PIPS mit der entgegennahme sehr strikt war und ich nichts in PIPS verändern wollte und keinem Umweg sah, ist das jetzt mit diesen globalen Variablen gelöst
    for(int i = 0; i <= nScenarios; i++) {
        myl.push_back(Linopy_Storage.at(0).Get_myl());//Diese beiden mussten außerhalb bestimmt werden. Weil sie nur in Block 0 zu finden sind.
        mzl.push_back(Linopy_Storage.at(0).Get_mzl());
    }
    for(int i = 0; i <= nScenarios; i++) {
        n.push_back(Linopy_Storage.at(i).Get_n());
        my.push_back(Linopy_Storage.at(i).Get_my());
        mz.push_back(Linopy_Storage.at(i).Get_mz());

        nnzA.push_back(Linopy_Storage.at(i).Get_nnz_A());
        nnzB.push_back(Linopy_Storage.at(i).Get_nnz_B());
        nnzC.push_back(Linopy_Storage.at(i).Get_nnz_C());
        nnzD.push_back(Linopy_Storage.at(i).Get_nnz_D());
        nnzBL.push_back(Linopy_Storage.at(i).Get_nnz_BL());
        nnzDL.push_back(Linopy_Storage.at(i).Get_nnz_DL());

        A_col.push_back(Linopy_Storage.at(i).Get_A_col());
        A_row.push_back(Linopy_Storage.at(i).Get_A_row());
        A_data.push_back(Linopy_Storage.at(i).Get_A_data());

        B_col.push_back(Linopy_Storage.at(i).Get_B_col());
        B_row.push_back(Linopy_Storage.at(i).Get_B_row());
        B_data.push_back(Linopy_Storage.at(i).Get_B_data());

        BL_col.push_back(Linopy_Storage.at(i).Get_BL_col());
        BL_row.push_back(Linopy_Storage.at(i).Get_BL_row());
        BL_data.push_back(Linopy_Storage.at(i).Get_BL_data());

        C_col.push_back(Linopy_Storage.at(i).Get_C_col());
        C_row.push_back(Linopy_Storage.at(i).Get_C_row());
        C_data.push_back(Linopy_Storage.at(i).Get_C_data());

        D_col.push_back(Linopy_Storage.at(i).Get_D_col());
        D_row.push_back(Linopy_Storage.at(i).Get_D_row());
        D_data.push_back(Linopy_Storage.at(i).Get_D_data());

        DL_col.push_back(Linopy_Storage.at(i).Get_DL_col());
        DL_row.push_back(Linopy_Storage.at(i).Get_DL_row());
        DL_data.push_back(Linopy_Storage.at(i).Get_DL_data());
    // Linopy stores the constraints and the objective function in vectors. 
        c.push_back(Linopy_Storage.at(i).Get_objective_function());  // Objective function constraints
        b.push_back(Linopy_Storage.at(i).Get_equality_vector());  // Equality Constraints
        bl.push_back(Linopy_Storage.at(i).Get_equality_linking_vector()); // Equality Linking Constraints

        dl.push_back(Linopy_Storage.at(i).Get_lower_inequality_constraint_vector()); // All kind of inequality Constraints
        idl.push_back(Linopy_Storage.at(i).Get_lower_inequality_constraint_vector_indicator()); //Indikator Vektor für dl. Nicht von linopy, für PIPS
        du.push_back(Linopy_Storage.at(i).Get_upper_inequality_constraint_vector());
        idu.push_back(Linopy_Storage.at(i).Get_upper_inequality_constraint_vector_indicator()); //Indikator Vektor für du. Nicht von linopy, für PIPS
        dll.push_back(Linopy_Storage.at(i).Get_lower_inequality_constraint_linking_vector());
        idll.push_back(Linopy_Storage.at(i).Get_lower_inequality_constraint_linking_vector_indicator()); //Indikator ...
        dlu.push_back(Linopy_Storage.at(i).Get_upper_inequality_constraint_linking_vector());
        idlu.push_back(Linopy_Storage.at(i).Get_upper_inequality_constraint_linking_vector_indicator()); //Indikator ..

        xl.push_back(Linopy_Storage.at(i).Get_lower_inequality_vector());
        ixl.push_back(Linopy_Storage.at(i).Get_lower_inequality_vector_indicator()); //Indikator ...
        xu.push_back(Linopy_Storage.at(i).Get_upper_inequality_vector());
        ixu.push_back(Linopy_Storage.at(i).Get_upper_inequality_vector_indicator()); //Indikator ...
    }
    
    // set callbacks
    FNNZ nCall = &nSize; // Number of Variables in a Block. For Linopy it would be x.
    FNNZ myCall = &mySize; // Number of equality constraints(no linking). Lenght of b.
    FNNZ mzCall = &mzSize; // Numer of inequality constraints(no linking). Lenght of dl(or du).
    FNNZ mylCall = &mylSize; // Number of equality constraints(linking). Lenght of bl.
    FNNZ mzlCall = &mzlSize; // Numer of inequality constraints(no linking). Lenght of dll(or dlu).
    FNNZ fnnzQ = &nnzAllZero; // Number of entries in Matrix Q, e.g. lenght of the data
    FNNZ fnnzA = &nnzMatEqStage1; // Number of entries in Matrix A, e.g. lenght of the data
    FNNZ fnnzB = &nnzMatEqStage2; // Number of entries in Matrix B, e.g. lenght of the data

    FNNZ fnnzC = &nnzMatIneqStage1; // Number of entries in Matrix C, e.g. lenght of the data
    FNNZ fnnzD = &nnzMatIneqStage2; // Number of entries in Matrix C, e.g. lenght of the data

    FNNZ fnnzBl = &nnzMatEqLink; // Number of entries in Matrix DL, e.g. lenght of the data
    FVEC fbl = &vecEqRhsLink; // Vector bl, equality linking constraint
    FMAT fBl = &matEqLink; // Matrix BL, equality Jacobian(Linking)
    FMAT fDl = &matIneqLink; // Matrix DL, inequality Jacobian(Linking)
    FNNZ fnnzDl = &nnzMatIneqLink; // nŃumber of entries in Matrix DL, e.g. lenght of the data

    FVEC fdlupp = &vecIneqRhsLink; // Vector dlu, upper inequality linking constraints
    FVEC fdllow = &vecIneqLhsLink; // Vector dll, upper inequality linking constraints
    FVEC fidlupp = &vecIneqRhsLinkActive; // Vector dlu indicator, if there is a constraint, e.g. no inf
    FVEC fidllow = &vecIneqLhsLinkActive; // Vector dll indicator, if there is a constraint, e.g. no inf

    FVEC fc = &vecObj; // Vector c, minimizing constraint
    FVEC fb = &vecEqRhs; // Vector b, equality constraint

    FVEC fclow = &vecIneqLhs; // Vector dl, lower inequality constraint
    FVEC fcupp = &vecIneqRhs;// Vector du, upper inequality constraint
    FVEC fxlow = &vecXlb; // Vector xl, lower value constraint
    FVEC fxupp = &vecXub; // Vector xu, upper value constraint
    FVEC ficlow = &vecIneqLhsActive; // Vector dl indicator, if there is a constraint, e.g. no inf
    FVEC fixlow = &vecXlbActive; // Vector xl indicator, if there is a constraint, e.g. no inf
    FVEC ficupp = &vecIneqRhsActive; // Vector du indicator, if there is a constraint, e.g. no inf
    FVEC fixupp = &vecXubActive; // Vector xu indicator, if there is a constraint, e.g. no inf


    FMAT fQ = &matAllZero; // Matrix Q, hessian
    FMAT fA = &matEqStage1; // Matrix A, equality Jacobian
    FMAT fB = &matEqStage2; // Matrix B, equality Jacobian
    FMAT fC = &matIneqStage1; // Matrix C, inequality Jacobian
    FMAT fD = &matIneqStage2; // Matrix D, inequality Jacobian
   
    ProbData probData(nScenarios);
    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //build the problem tree
    std::unique_ptr<DistributedInputTree::DistributedInputNode> data_root = std::make_unique<DistributedInputTree::DistributedInputNode>(&probData, 0, nCall, myCall, mylCall, mzCall, mzlCall, fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB, fBl,
         fnnzBl, fb, fbl, fC, fnnzC, fD, fnnzD, fDl, fnnzDl, fclow, ficlow, fcupp, ficupp, fdllow, fidllow, fdlupp, fidlupp, fxlow, fixlow, fxupp,
         fixupp, nullptr, false);

    auto* root = new DistributedInputTree(std::move(data_root));
    for (int id = 1; id <= nScenarios; id++) {
        std::unique_ptr<DistributedInputTree::DistributedInputNode> data_child = std::make_unique<DistributedInputTree::DistributedInputNode>(&probData, id, nCall, myCall, mylCall, mzCall, mzlCall, fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB,
            fBl, fnnzBl, fb, fbl, fC, fnnzC, fD, fnnzD, fDl, fnnzDl, fclow, ficlow, fcupp, ficupp, fdllow, fidllow, fdlupp, fidlupp, fxlow, fixlow,
            fxupp, fixupp, nullptr, false);
      
        root->add_child(std::make_unique<DistributedInputTree>(std::move(data_child)));
    }

    if (rank == 0)
        std::cout << "Using a total of " << size << " MPI processes.\n";
   
    /* use BiCGStab for outer solve */
    pipsipmpp_options::set_int_parameter("INNER_SC_SOLVE", 0);
    std::cout << "PIPSIPMppInterface pipsIpm" << std::endl;

    PIPSIPMppInterface pipsIpm(root, InteriorPointMethodType::PRIMAL, MPI_COMM_WORLD, ScalerType::GEOMETRIC_MEAN, PresolverType::PRESOLVE); 

    if (rank == 0)
        std::cout << "PIPSIPMppInterface created\n";

    if (rank == 0)
        std::cout << "solving...\n";


    pipsIpm.run();

    const double objective = pipsIpm.getObjective();
    if (rank == 0)
        std::cout << "solving finished ... objective value: " << objective << "\n";

    delete root;

    MPI_Finalize();

    return 0;
}
