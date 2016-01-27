//
// boilerplate.h
//
// Common code for creating Wolfram LibraryLink DLLs
//
#pragma once

#include "windows.h"
#include "WolframLibrary.h"
#include "memory.h"
#include "strsafe.h"

void Trace(LPSTR format, ...)
    {
    char message[256];
    va_list args = NULL;
    va_start(args, format);
    StringCchVPrintfA(message, 256, format, args);

    char buffer[256];
    StringCchCopyA(buffer, 256, "ode: ");
    StringCchCatA(buffer, 256, message);
    OutputDebugStringA(buffer);
    }

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

#include "odeStep.h"

ODESTATE runODECore(const ODESTATE& inits, const ODEPARAMETERS & params, double tDuration, double dt)
    {
    double t = 0;
    ODESTATE stateA = inits;
    ODESTATE stateB;
    ODESTATE* pStateCur  = &stateA;
    ODESTATE* pStateNext = &stateB;

    if (true)
        {
        while (t <= tDuration)
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
    
mint arraySize(WolframLibraryData libData, MTensor tensor)
    {
    mint const* dims = libData->MTensor_getDimensions(tensor);
    return dims[0];
    }
void copyFromTensor(WolframLibraryData libData, double* pDest, mint cDest, MTensor tensor)
    {
    double* pSource = libData->MTensor_getRealData(tensor);

    mint cSrc = arraySize(libData, tensor);
    mint cb = min(cDest, cSrc) * sizeof(double);
    mint cbDest = cDest * sizeof(double);

    memcpy(pDest, pSource, cb);
    memset(&pDest[cSrc], 0, cbDest - cb);
    }
void copyToTensor(WolframLibraryData libData, MTensor tensor, double* pSource, mint cSrc)
    {
    double* pDest = libData->MTensor_getRealData(tensor);

    mint cDest = arraySize(libData, tensor);
    mint cb = min(cDest, cSrc) * sizeof(double);
    mint cbDest = cDest * sizeof(double);

    memcpy(pDest, pSource, cb);
    memset(&pDest[cSrc], 0, cbDest - cb);
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
    double  dt           = MArgument_getReal(argv[3]);

    mint rankState  = (mint)(sizeof(ODESTATE)      / sizeof(double));
    mint rankParams = (mint)(sizeof(ODEPARAMETERS) / sizeof(double));

    Trace("state:%d params:%d", rankState, rankParams);

    ODESTATE inits;
    ODEPARAMETERS parameters;

    copyFromTensor(libData, (double*)&inits, rankState, tensorInits);
    copyFromTensor(libData, (double*)&parameters, rankParams, tensorParams);
    
    ODESTATE result = runODECore(inits, parameters, duration, dt);

    MTensor tensorResult;
    libData->MTensor_new(MType_Real, 1, &rankState, &tensorResult);
    copyToTensor(libData, tensorResult, (double*)&result, rankState);

    MArgument_setMTensor(aResult, tensorResult);

    return LIBRARY_NO_ERROR;
    }