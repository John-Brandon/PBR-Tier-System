#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=/usr/local/bin/gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/BRENT.o \
	${OBJECTDIR}/Declare_variables_module.o \
	${OBJECTDIR}/Format_module.o \
	${OBJECTDIR}/Generate_random_numbers_module.o \
	${OBJECTDIR}/Initialize_pop_module.o \
	${OBJECTDIR}/PBR_Errorcheck_module.o \
	${OBJECTDIR}/PBR_FileIO_Module.o \
	${OBJECTDIR}/PBRmodule.o \
	${OBJECTDIR}/Random_module.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pbr_netbeans

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pbr_netbeans: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pbr_netbeans ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/BRENT.o: BRENT.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/BRENT.o BRENT.f90

${OBJECTDIR}/Declare_variables_module.o: Declare_variables_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/Declare_variables_module.o Declare_variables_module.f90

${OBJECTDIR}/Format_module.o: Format_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/Format_module.o Format_module.f90

${OBJECTDIR}/Generate_random_numbers_module.o: Generate_random_numbers_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/Generate_random_numbers_module.o Generate_random_numbers_module.f90

${OBJECTDIR}/Initialize_pop_module.o: Initialize_pop_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/Initialize_pop_module.o Initialize_pop_module.f90

${OBJECTDIR}/PBR_Errorcheck_module.o: PBR_Errorcheck_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/PBR_Errorcheck_module.o PBR_Errorcheck_module.f90

${OBJECTDIR}/PBR_FileIO_Module.o: PBR_FileIO_Module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/PBR_FileIO_Module.o PBR_FileIO_Module.f90

${OBJECTDIR}/PBRmodule.o: PBRmodule.f 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/PBRmodule.o PBRmodule.f

${OBJECTDIR}/Random_module.o: Random_module.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/Random_module.o Random_module.f90

${OBJECTDIR}/main.o: main.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/main.o main.f90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/pbr_netbeans
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
