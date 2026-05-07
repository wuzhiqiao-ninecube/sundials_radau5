/* ---------------------------------------------------------------------------
 * python_stubs.c — Stubs for SUNDIALS Python FunctionTable_Destroy symbols
 *
 * When SUNDIALS is built with SUNDIALS_ENABLE_PYTHON, the core library
 * references these symbols which are normally provided by the Python
 * bindings layer.  For standalone C programs we provide no-op stubs.
 * ---------------------------------------------------------------------------*/

void SUNContextFunctionTable_Destroy(void* ptr) { (void)ptr; }
void SUNLinearSolverFunctionTable_Destroy(void* ptr) { (void)ptr; }
void SUNNonlinearSolverFunctionTable_Destroy(void* ptr) { (void)ptr; }
void SUNDomEigEstimatorFunctionTable_Destroy(void* ptr) { (void)ptr; }
void SUNStepperFunctionTable_Destroy(void* ptr) { (void)ptr; }
