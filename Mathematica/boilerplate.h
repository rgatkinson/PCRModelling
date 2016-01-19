//
// boilerplate.h
//
// Common code for creating Wolfram LibraryLink DLLs
//
#pragma once

#include "WolframLibrary.h"
#include "memory.h"

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion() 
    {
    return WolframLibraryVersion;
    }

EXTERN_C DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData)
    {
    return 0;
    }

EXTERN_C DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) 
    {
    return;
    }

EXTERN_C DLLEXPORT int demo_I_I(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) 
    {
    mint I0;
    mint I1;
    I0 = MArgument_getInteger(Args[0]);
    I1 = I0 + 2;
    MArgument_setInteger(Res, I1);
    return LIBRARY_NO_ERROR;
    }

#include "odeStep.h"

ODESTATE runODECore(const ODESTATE& inits, const ODEPARAMETERS & params, double tDuration)
    {
    const int nSteps = 1000;
    const double dt = tDuration / nSteps;

    double t = 0;
    ODESTATE stateA = inits;
    ODESTATE stateB;
    ODESTATE* pStateCur  = &stateA;
    ODESTATE* pStateNext = &stateB;

    if (false)
        {
        for (int iStep = 0; iStep < nSteps; iStep++)
            {
            // Run one step of the ODE
            odeStep(t, dt, params, *pStateCur, *pStateNext);

            // Advance time
            t += dt;

            // Swap state locations for the next iteration
            ODESTATE* pT = pStateCur;
            pStateCur = pStateNext;
            pStateNext = pT;
            }
        }

    // Result is in pStateCur
    return *pStateCur;
    }
    
void copyFromTensor(WolframLibraryData libData, mint count, double* pDest, MTensor tensor)
    {
    double* pSource = libData->MTensor_getRealData(tensor);
    memcpy(pDest, pSource, count * sizeof(double));
    }
void copyToTensor(WolframLibraryData libData, mint count, MTensor tensor, double* pSource)
    {
    double* pDest = libData->MTensor_getRealData(tensor);
    memcpy(pDest, pSource, count * sizeof(double));
    }

// Parameter order:
// inits
// params
// duration
EXTERN_C DLLEXPORT int runODE(WolframLibraryData libData, mint argc, MArgument *argv, MArgument aResult)
    {
    MTensor tensorInits  = MArgument_getMTensor(argv[0]);
    MTensor tensorParams = MArgument_getMTensor(argv[1]);
    double  duration     = MArgument_getReal(argv[2]);

    mint rankState  = (mint)(sizeof(ODESTATE)      / sizeof(double));
    mint rankParams = (mint)(sizeof(ODEPARAMETERS) / sizeof(double));

    ODESTATE inits;
    ODEPARAMETERS parameters;

    copyFromTensor(libData, rankState, (double*)&inits, tensorInits);
    copyFromTensor(libData, rankParams, (double*)&parameters, tensorParams);
    
    ODESTATE result = runODECore(inits, parameters, duration);

    MTensor tensorResult;
    libData->MTensor_new(MType_Real, 1, &rankState, &tensorResult);
    copyToTensor(libData, rankState, tensorResult, (double*)&result);

    MArgument_setMTensor(aResult, tensorResult);

    return LIBRARY_NO_ERROR;
    }