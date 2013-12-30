EXE                    = ../bin/nomad
EXE_MPI		           = ../bin/nomad.MPI
LIB                    = ../lib/nomad.a
LIB_MPI                = ../lib/nomad.MPI.a
COMPILATOR             = g++
COMPILATOR_MPI = mpic++
SUNAME = $(shell uname)
OSS=$(findstring MINGW32,$(SUNAME)) 
ifneq "$(strip $(OSS))" ""
COMPILATOR_MPI = g++
EXE                    = ../bin/nomad.exe
EXE_MPI		           = ../bin/nomad.MPI.exe
endif

BUILD_DIR_NO_MPI = ../builds/ObjsNoMPI 
BUILD_DIR_MPI = ../builds/ObjsMPI
BUILD_DIR = 
COMPILATOR_OPTIONS     = -O2 -ansi 
# COMPILATOR_OPTIONS = -g -ansi -Wall 
# COMPILATOR_OPTIONS_MPI  = $(COMPILATOR_OPTIONS) -DUSE_MPI
COMPILATOR_OPTIONS_MPI = $(COMPILATOR_OPTIONS) -DUSE_MPI # -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
LIBS                   = -lm
LIBS_MPI               = $(LIBS) -lmpi
INCLUDE                = -I.
COMPILE                = 
COMPILE_NO_MPI         = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
COMPILE_MPI            = $(COMPILATOR_MPI) $(COMPILATOR_OPTIONS_MPI) $(INCLUDE) -c
OBJS_LIB               = $(BUILD_DIR)/Barrier.o $(BUILD_DIR)/Cache.o $(BUILD_DIR)/Cache_File_Point.o $(BUILD_DIR)/Cache_Point.o \
                         $(BUILD_DIR)/Cache_Search.o $(BUILD_DIR)/Clock.o $(BUILD_DIR)/Direction.o $(BUILD_DIR)/Directions.o $(BUILD_DIR)/Display.o \
                         $(BUILD_DIR)/Double.o $(BUILD_DIR)/Eval_Point.o $(BUILD_DIR)/Evaluator.o $(BUILD_DIR)/Evaluator_Control.o \
                         $(BUILD_DIR)/Exception.o $(BUILD_DIR)/Extended_Poll.o $(BUILD_DIR)/L_Curve.o $(BUILD_DIR)/LH_Search.o $(BUILD_DIR)/Mads.o \
                         $(BUILD_DIR)/Mesh.o $(BUILD_DIR)/Model_Sorted_Point.o $(BUILD_DIR)/Model_Stats.o $(BUILD_DIR)/Multi_Obj_Evaluator.o \
                         $(BUILD_DIR)/Parameters.o $(BUILD_DIR)/Parameter_Entries.o $(BUILD_DIR)/Parameter_Entry.o \
                         $(BUILD_DIR)/Pareto_Front.o $(BUILD_DIR)/Pareto_Point.o $(BUILD_DIR)/Phase_One_Evaluator.o \
                         $(BUILD_DIR)/Phase_One_Search.o $(BUILD_DIR)/Point.o $(BUILD_DIR)/Priority_Eval_Point.o $(BUILD_DIR)/Quad_Model.o \
                         $(BUILD_DIR)/Quad_Model_Evaluator.o $(BUILD_DIR)/Quad_Model_Search.o $(BUILD_DIR)/Random_Pickup.o \
                         $(BUILD_DIR)/RNG.o $(BUILD_DIR)/Signature.o $(BUILD_DIR)/Slave.o $(BUILD_DIR)/Speculative_Search.o $(BUILD_DIR)/Stats.o \
                         $(BUILD_DIR)/utils.o $(BUILD_DIR)/Variable_Group.o $(BUILD_DIR)/VNS_Search.o
                         
OBJS                   = $(BUILD_DIR)/nomad.o $(OBJS_LIB)


all:
	@mkdir -p $(BUILD_DIR_NO_MPI) 
	$(MAKE) BUILD_DIR=$(BUILD_DIR_NO_MPI) COMPILE='$(COMPILE_NO_MPI)' $(EXE) $(LIB)
	
mpi:
	@mkdir -p $(BUILD_DIR_MPI)
	$(MAKE) BUILD_DIR=$(BUILD_DIR_MPI) COMPILE='$(COMPILE_MPI)' $(EXE_MPI) $(LIB_MPI)

installAllNomad: all mpi
	@echo "cleaning obj files and build directories"
	@rm -rf $(BUILD_DIR_NO_MPI) $(BUILD_DIR_MPI)
	
clean:
	@echo "cleaning obj files and build directories"
	@rm -rf $(BUILD_DIR_NO_MPI) $(BUILD_DIR_MPI)
	
del: clean
	@echo "cleaning trash files"
	@rm -f core *~
	@echo "cleaning exe file"
	@rm -f $(EXE) $(EXE_MPI)
	@echo "cleaning lib files"
	@rm -f $(LIB) $(LIB_MPI)	
	
$(EXE): $(OBJS)
	@echo "building the release version (no MPI)"
	@echo "exe file: "$(EXE)
	@$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)
	@strip $(EXE)

$(LIB): $(OBJS_LIB)	
	@echo "building the library (no MPI)"
	@echo "lib file: "$(LIB)
	@rm -f $(LIB)
	@ar -r $(LIB) $(OBJS_LIB)
	@ranlib $(LIB)
 
$(EXE_MPI): $(OBJS)
	@echo "building the release version (MPI)"
	@echo "exe file: "$(EXE_MPI)
	@$(COMPILATOR_MPI) -o $(EXE_MPI) $(OBJS) $(LIBS_MPI) $(COMPILATOR_OPTIONS_MPI)	
	@strip $(EXE_MPI)

$(LIB_MPI): $(OBJS_LIB)
	@echo "building the MPI library"
	@echo "lib file: "$(LIB_MPI)
	@rm -f $(LIB_MPI)
	@ar -r $(LIB_MPI) $(OBJS_LIB)
	@ranlib $(LIB_MPI)
	
#$(BUILD_DIR)/%.o: %.cpp %.hpp
#	$(COMPILE) $< -o $@
	
$(BUILD_DIR)/Barrier.o: Barrier.cpp Barrier.hpp Filter_Point.hpp Set_Element.hpp
	$(COMPILE) Barrier.cpp -o $@

$(BUILD_DIR)/Cache.o: Cache.cpp Cache.hpp
	$(COMPILE) Cache.cpp -o $@

$(BUILD_DIR)/Cache_File_Point.o: Cache_File_Point.cpp Cache_File_Point.hpp \
                    Eval_Point.cpp Eval_Point.hpp Uncopyable.hpp
	$(COMPILE) Cache_File_Point.cpp -o $@

$(BUILD_DIR)/Cache_Point.o: Cache_Point.cpp Cache_Point.hpp
	$(COMPILE) Cache_Point.cpp -o $@

$(BUILD_DIR)/Cache_Search.o: Cache_Search.cpp Cache_Search.hpp Search.hpp
	$(COMPILE) Cache_Search.cpp -o $@

$(BUILD_DIR)/Clock.o: Clock.cpp Clock.hpp
	$(COMPILE) Clock.cpp -o $@

$(BUILD_DIR)/Direction.o: Direction.cpp Direction.hpp Point.hpp
	$(COMPILE) Direction.cpp -o $@

$(BUILD_DIR)/Directions.o: Directions.cpp Directions.hpp Mesh.hpp Random_Pickup.hpp RNG.hpp
	$(COMPILE) Directions.cpp -o $@

$(BUILD_DIR)/Display.o: Display.cpp Display.hpp utils.hpp
	$(COMPILE) Display.cpp -o $@

$(BUILD_DIR)/Double.o: Double.cpp Double.hpp Exception.hpp Display.hpp
	$(COMPILE) Double.cpp -o $@

$(BUILD_DIR)/Eval_Point.o: Eval_Point.cpp Eval_Point.hpp Parameters.hpp Cache_File_Point.hpp \
              Set_Element.hpp
	$(COMPILE) Eval_Point.cpp -o $@

$(BUILD_DIR)/Evaluator.o: Evaluator.cpp Evaluator.hpp Priority_Eval_Point.hpp Stats.hpp
	$(COMPILE) Evaluator.cpp -o $@

$(BUILD_DIR)/Evaluator_Control.o: Evaluator_Control.cpp Evaluator_Control.hpp \
            Barrier.hpp Pareto_Front.hpp Slave.hpp Quad_Model.hpp TGP_Model.hpp        
	$(COMPILE) Evaluator_Control.cpp -o $@

$(BUILD_DIR)/Exception.o: Exception.cpp Exception.hpp
	$(COMPILE) Exception.cpp -o $@

$(BUILD_DIR)/Extended_Poll.o: Extended_Poll.cpp Extended_Poll.hpp Signature_Element.hpp \
                 Set_Element.hpp Mads.hpp
	$(COMPILE) Extended_Poll.cpp -o $@

$(BUILD_DIR)/L_Curve.o: L_Curve.cpp L_Curve.hpp Double.hpp Uncopyable.hpp
	$(COMPILE) L_Curve.cpp -o $@

$(BUILD_DIR)/LH_Search.o: LH_Search.cpp LH_Search.hpp Search.hpp Mads.hpp RNG.hpp Evaluator_Control.hpp
	$(COMPILE) LH_Search.cpp -o $@

$(BUILD_DIR)/Mads.o: Mads.cpp Mads.hpp Evaluator_Control.hpp L_Curve.hpp \
        LH_Search.hpp LH_Search.cpp \
        Speculative_Search.cpp Speculative_Search.hpp \
        Extended_Poll.cpp Extended_Poll.hpp \
        VNS_Search.hpp VNS_Search.cpp \
        Quad_Model_Search.hpp Quad_Model_Search.cpp \
        Cache_Search.hpp Cache_Search.cpp \
        Phase_One_Search.cpp Phase_One_Search.hpp
	$(COMPILE) Mads.cpp -o $@

$(BUILD_DIR)/Mesh.o: Mesh.cpp Mesh.hpp Point.hpp
	$(COMPILE) Mesh.cpp -o $@

$(BUILD_DIR)/Multi_Obj_Evaluator.o: Multi_Obj_Evaluator.cpp Multi_Obj_Evaluator.hpp Phase_One_Evaluator.hpp
	$(COMPILE) Multi_Obj_Evaluator.cpp -o $@

$(BUILD_DIR)/Model_Sorted_Point.o: Model_Sorted_Point.cpp Model_Sorted_Point.hpp Point.hpp
	$(COMPILE) Model_Sorted_Point.cpp -o $@

$(BUILD_DIR)/Model_Stats.o: Model_Stats.cpp Model_Stats.hpp Double.hpp
	$(COMPILE) Model_Stats.cpp -o $@

$(BUILD_DIR)/nomad.o: nomad.cpp nomad.hpp  Mads.hpp
	$(COMPILE) nomad.cpp -o $@

$(BUILD_DIR)/Parameters.o: Parameters.cpp Parameters.hpp  Parameter_Entries.hpp Signature.hpp
	$(COMPILE) Parameters.cpp -o $@

$(BUILD_DIR)/Parameter_Entries.o: Parameter_Entries.cpp Parameter_Entries.hpp Parameter_Entry.hpp
	$(COMPILE) Parameter_Entries.cpp -o $@

$(BUILD_DIR)/Parameter_Entry.o: Parameter_Entry.hpp Parameter_Entry.cpp  Display.hpp Uncopyable.hpp
	$(COMPILE) Parameter_Entry.cpp -o $@

$(BUILD_DIR)/Pareto_Front.o: Pareto_Front.cpp Pareto_Front.hpp Pareto_Point.hpp
	$(COMPILE) Pareto_Front.cpp -o $@

$(BUILD_DIR)/Pareto_Point.o: Pareto_Point.cpp Pareto_Point.hpp Multi_Obj_Evaluator.hpp
	$(COMPILE) Pareto_Point.cpp -o $@

$(BUILD_DIR)/Phase_One_Evaluator.o: Phase_One_Evaluator.cpp Phase_One_Evaluator.hpp Evaluator.hpp
	$(COMPILE) Phase_One_Evaluator.cpp -o $@

$(BUILD_DIR)/Phase_One_Search.o: Phase_One_Search.cpp Phase_One_Search.hpp Mads.hpp \
                    Search.hpp Evaluator_Control.hpp
	$(COMPILE) Phase_One_Search.cpp -o $@

$(BUILD_DIR)/Point.o: Point.cpp Point.hpp Double.hpp
	$(COMPILE) Point.cpp -o $@

$(BUILD_DIR)/Priority_Eval_Point.o: Priority_Eval_Point.cpp Priority_Eval_Point.hpp Eval_Point.hpp \
                       Set_Element.hpp
	$(COMPILE) Priority_Eval_Point.cpp -o $@

$(BUILD_DIR)/Quad_Model.o: Quad_Model.cpp Quad_Model.hpp Cache.hpp Model_Sorted_Point.hpp
	$(COMPILE) Quad_Model.cpp -o $@

$(BUILD_DIR)/Quad_Model_Evaluator.o: Quad_Model_Evaluator.cpp Quad_Model_Evaluator.hpp \
                        Evaluator.hpp Search.hpp
	$(COMPILE) Quad_Model_Evaluator.cpp -o $@

$(BUILD_DIR)/Quad_Model_Search.o: Quad_Model_Search.cpp Quad_Model_Search.hpp Mads.hpp \
                     Quad_Model_Evaluator.hpp
	$(COMPILE) Quad_Model_Search.cpp -o $@

$(BUILD_DIR)/Random_Pickup.o: Random_Pickup.cpp Random_Pickup.hpp RNG.cpp RNG.hpp Uncopyable.hpp
	$(COMPILE) Random_Pickup.cpp -o $@

$(BUILD_DIR)/RNG.o: RNG.cpp RNG.hpp defines.hpp
	$(COMPILE) RNG.cpp -o $@

$(BUILD_DIR)/Signature.o: Signature.cpp Signature.hpp Variable_Group.hpp
	$(COMPILE) Signature.cpp -o $@

$(BUILD_DIR)/Slave.o: Slave.cpp Slave.hpp Evaluator.hpp
	$(COMPILE) Slave.cpp -o $@

$(BUILD_DIR)/Speculative_Search.o: Speculative_Search.cpp Speculative_Search.hpp Mads.hpp Search.hpp \
                      Evaluator_Control.hpp
	$(COMPILE) Speculative_Search.cpp -o $@

$(BUILD_DIR)/Stats.o: Stats.cpp Stats.hpp  Clock.hpp Double.hpp Model_Stats.hpp
	$(COMPILE) Stats.cpp -o $@

TGP_Model.o: TGP_Model.cpp TGP_Model.hpp TGP_Output_Model.hpp Cache.hpp Model_Sorted_Point.hpp
	$(COMPILE) TGP_Model.cpp

TGP_Model_Evaluator.o: TGP_Model_Evaluator.cpp TGP_Model_Evaluator.hpp \
                       Evaluator.hpp Search.hpp
	$(COMPILE) TGP_Model_Evaluator.cpp

TGP_Model_Search.o: TGP_Model_Search.cpp TGP_Model_Search.hpp \
                    LH_Search.hpp TGP_Model_Evaluator.hpp
	$(COMPILE) TGP_Model_Search.cpp

TGP_Output_Model.o: TGP_Output_Model.cpp TGP_Output_Model.hpp Eval_Point.o $(LIB_TGP)
	$(COMPILE) TGP_Output_Model.cpp

$(BUILD_DIR)/utils.o: utils.cpp utils.hpp defines.hpp 
	$(COMPILE) utils.cpp -o $@

$(BUILD_DIR)/Variable_Group.o: Variable_Group.cpp Variable_Group.hpp  Directions.hpp
	$(COMPILE) Variable_Group.cpp -o $@

$(BUILD_DIR)/VNS_Search.o: VNS_Search.cpp VNS_Search.hpp Search.hpp Evaluator_Control.hpp
	$(COMPILE) VNS_Search.cpp -o $@


