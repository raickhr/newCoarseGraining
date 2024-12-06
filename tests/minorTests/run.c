if (solve_normal) {
     PetscCall(KSPSetOperators(ksp, N, N));
   } 
   PetscCall(KSPSetInitialGuessNonzero(ksp, nonzero_guess));
   PetscCall(KSPSetFromOptions(ksp));

   /*
      Here we explicitly call KSPSetUp() and KSPSetUpOnBlocks() to
      enable more precise profiling of setting up the preconditioner.
      These calls are optional, since both will be called within
      KSPSolve() if they haven't been called already.
   */
   PetscCall(KSPSetUp(ksp));
   PetscCall(KSPSetUpOnBlocks(ksp));

   /*
                          Solve system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   /*
      Begin profiling next stage
   */
   PetscPreLoadStage("KSPSolve");

   /*
      Solve linear system
   */
   if (solve_normal) {
     PetscCall(KSPSolve(ksp, Ab, x));
   }
   PetscCall(PetscObjectSetName((PetscObject)x, "x"));

   /*
       Conclude profiling this stage
    