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
    


    std::vector<int*> A_col;
    std::vector<int*> A_row;
    std::vector<double*> A_data;

    std::vector<int*> B_col;
    std::vector<int*> B_row;
    std::vector<double*> B_data;

    std::vector<int*> BL_col;
    std::vector<int*> BL_row;
    std::vector<double*> BL_data;

    std::vector<int*> C_col;
    std::vector<int*> C_row;
    std::vector<double*> C_data;

    std::vector<int*> D_col;
    std::vector<int*> D_row;
    std::vector<double*> D_data;

    std::vector<int*> DL_col;
    std::vector<int*> DL_row;
    std::vector<double*> DL_data;

    //Linopy stores the constraints and the objective function in vectors. 
    std::vector<double*> c;  // Objective function constraints
    std::vector<double*> b;  // Equality Constraints
    std::vector<double*> bl; // Equality Linking Constraints

    std::vector<double*> dl; // All kind of inequality Constraints
    std::vector<double*> idl; //Indikator Vektor für dl. Nicht von linopy, für PIPS
    std::vector<double*> du;
    std::vector<double*> idu; //Indikator Vektor für du. Nicht von linopy, für PIPS
    std::vector<double*> dll;
    std::vector<double*> idll; //Indikator ...
    std::vector<double*> dlu;
    std::vector<double*> idlu; //Indikator ..

    std::vector<double*> xl;
    std::vector<double*> ixl; //Indikator ...
    std::vector<double*> xu;
    std::vector<double*> ixu; //Indikator ...


    std::vector<int> n;
    std::vector<int> my;
    std::vector<int> myl;
    std::vector<int> mz;
    std::vector<int> mzl;

    std::vector<int> nnzA;
    std::vector<int> nnzB;
    std::vector<int> nnzC;
    std::vector<int> nnzD;
    std::vector<int> nnzBL;
    std::vector<int> nnzDL;
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
        *nnz = myl.at(id);
        return 0;
    }
    int mzlSize(void*, int id, int* nnz) {
        *nnz = mzl.at(id);
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
        vec = std::move(b.at(id));
        return 0;
    }
    int vecEqRhsLink(void*, int id, double* vec, int len) {
        vec = std::move(bl.at(0));
        return 0;
    }
    int vecIneqRhs(void*, int id, double* vec, int len) {
        vec = std::move(du.at(id));
        return 0;
    }
    int vecIneqRhsActive(void*, int id, double* vec, int len) {
        vec = std::move(idu.at(id));
        return 0;
    }
    int vecIneqLhs(void*, int id, double* vec, int len) {
        vec = std::move(dl.at(id));
        return 0;
    }
    int vecIneqLhsActive(void*, int id, double* vec, int len) {
        vec = std::move(idl.at(id));
        return 0;
    }
    int vecIneqRhsLink(void*, int id, double* vec, int len) {
        vec = std::move(dlu.at(0));
        return 0;
    }
    int vecIneqRhsLinkActive(void*, int id, double* vec, int len) {
        vec = std::move(idlu.at(0));
        return 0;
    }
    int vecIneqLhsLink(void*, int id, double* vec, int len) {
        vec = std::move(dll.at(id));
        return 0;
    }
    int vecIneqLhsLinkActive(void*, int id, double* vec, int len) {
        vec = std::move(idll.at(id));
        return 0;
    }
    int vecObj(void*, int id, double* vec, int len) {
        vec = std::move(c.at(id));
        return 0;
    }
    int vecXlb(void*, int id, double* vec, int len) {
        vec = std::move(xl.at(id));
        return 0;
    }
    int vecXub(void*, int id, double* vec, int len) {
        vec = std::move(xu.at(id));
        return 0;
    }
    int vecXlbActive(void*, int id, double* vec, int len) {
        vec = std::move(ixl.at(id));
        return 0;
    }
    int vecXubActive(void*, int id, double* vec, int len) {
        vec = std::move(ixu.at(id));
        return 0;
    }
    int matAllZero(void*, int, int* krowM, int* , double*) {
        krowM[0] = 0;
        krowM[1] = 0;
        return 0;
    }
    int matEqStage1(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(A_row.at(id));
        jcolM = std::move(A_col.at(id));
        M = std::move(A_data.at(id));
        return 0;
    }
    int matIneqLink(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(DL_row.at(id));
        jcolM = std::move(DL_col.at(id));
        M = std::move(DL_data.at(id));
        return 0;
    }
    int matIneqStage1(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(C_row.at(id));
        jcolM = std::move(C_col.at(id));
        M = std::move(C_data.at(id));
        std::cout << "C " <<krowM[5] << std::endl;
        return 0;
    }
    int matEqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(B_row.at(id));
        jcolM = std::move(B_col.at(id));
        M = std::move(B_data.at(id));
        return 0;
    }
    int matIneqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(D_row.at(id));
        jcolM = std::move(D_col.at(id));
        M = std::move(D_data.at(id));
        return 0;
    }
    int matEqLink(void*, int id, int* krowM, int* jcolM, double* M) {
        krowM = std::move(BL_row.at(id));
        jcolM = std::move(BL_col.at(id));
        M = std::move(BL_data.at(id));
        return 0;
    }


int main(int argc, char** argv) {
   
    // This is for the Parallel Processor
    MPI_Init(&argc, &argv);
   
    // This has to be the same number as the number of Blocks you put in Linopy.
    // Linopy gives an extra block out. Substract 1 from it.
    int nScenarios = 5; //Argument aus der Main-Funktion

    // All of these are function pointers. The type of function is defined above.
    // FNNZ is for a function which returns an integer to PIPS
    // FVEC is a function which returns an Vector to PIPS
    // FMAT is a function which returns a Matric to PIPS
    // All of those values are returned through pointers in the function arguments


   //std::cout<< "kommen wir hierhin" << std::endl;
    //Linopy_Init Model("/home/ken/Desktop/pypsa_with_PIPS/pypsa-model",nScenarios);
    //std::cout << "Was passiert hier" << std::endl;

    std::vector<Linopy> Linopy_Storage;
    for(int i = 0; i <= nScenarios; i++) {
        Linopy_Storage.push_back(Linopy("/home/ken/Desktop/pypsa_with_PIPS/pypsa-model", i, nScenarios));
        std::cout << "Linopy_Storage " << i << " wurde kreirt" << std::endl;
    }

    Linopy_Storage.at(0).Transform_Matrix_Cols();
    std::cout << "Hier wurde Transform_Matrix_Cols ausgeführt" << std::endl;
    for(int i = 1; i <= nScenarios; i++) {
        Linopy_Storage.at(i).Set_xvec_A(Linopy_Storage.at(0).Get_xvec_A(),Linopy_Storage.at(0).Get_xvec_A_size());
        Linopy_Storage.at(i).Transform_Matrix_Cols();
    }

    for(int i = 0; i <= nScenarios; i++) {
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

        n.push_back(Linopy_Storage.at(i).Get_n());
        my.push_back(Linopy_Storage.at(i).Get_my());
        myl.push_back(Linopy_Storage.at(i).Get_myl());
        mz.push_back(Linopy_Storage.at(i).Get_mz());
        mzl.push_back(Linopy_Storage.at(i).Get_mzl());

        nnzA.push_back(Linopy_Storage.at(i).Get_nnz_A());
        nnzB.push_back(Linopy_Storage.at(i).Get_nnz_B());
        nnzC.push_back(Linopy_Storage.at(i).Get_nnz_C());
        nnzD.push_back(Linopy_Storage.at(i).Get_nnz_D());
        nnzBL.push_back(Linopy_Storage.at(i).Get_nnz_BL());
        nnzDL.push_back(Linopy_Storage.at(i).Get_nnz_DL());
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
