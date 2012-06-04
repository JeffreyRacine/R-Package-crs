EXE                    = ../bin/nomad.exe
EXE_MPI		       = ../bin/nomad.MPI.exe
LIB                    = ../lib/nomad.a
LIB_MPI                = ../lib/nomad.MPI.a
COMPILATOR             = g++
# COMPILATOR_NO_MPI    = g++
COMPILATOR_MPI = mpic++
SUNAME = $(shell uname)
OSS=$(findstring MINGW32,$(SUNAME)) 
ifneq "$(strip $(OSS))" ""
COMPILATOR_MPI		= g++
endif

BUILD_DIR_NO_MPI = ObjsNoMPI 
BUILD_DIR_MPI = ObjsMPI
BUILD_DIR = 
COMPILATOR_OPTIONS     = -O2 -ansi
COMPILATOR_OPTIONS_MPI = $(COMPILATOR_OPTIONS) -DUSE_MPI
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
	
$(BUILD_DIR)/Barrier.o: Barrier.cpp Barrier.hpp Filter_Point.hpp Set_Element.hpp $(BUILD_DIR)/Eval_Point.o $(BUILD_DIR)/Cache.o
	$(COMPILE) Barrier.cpp -o $@

$(BUILD_DIR)/Cache.o: Cache.cpp Cache.hpp $(BUILD_DIR)/Cache_Point.o $(BUILD_DIR)/Clock.o
	$(COMPILE) Cache.cpp -o $@

$(BUILD_DIR)/Cache_File_Point.o: Cache_File_Point.cpp Cache_File_Point.hpp $(BUILD_DIR)/Display.o \
                    Eval_Point.cpp Eval_Point.hpp Uncopyable.hpp
	$(COMPILE) Cache_File_Point.cpp -o $@

$(BUILD_DIR)/Cache_Point.o: Cache_Point.cpp Cache_Point.hpp $(BUILD_DIR)/Eval_Point.o
	$(COMPILE) Cache_Point.cpp -o $@

$(BUILD_DIR)/Cache_Search.o: Cache_Search.cpp Cache_Search.hpp $(BUILD_DIR)/Mads.o Search.hpp $(BUILD_DIR)/Evaluator_Control.o
	$(COMPILE) Cache_Search.cpp -o $@

$(BUILD_DIR)/Clock.o: Clock.cpp Clock.hpp
	$(COMPILE) Clock.cpp -o $@

$(BUILD_DIR)/Direction.o: Direction.cpp Direction.hpp $(BUILD_DIR)/Point.o
	$(COMPILE) Direction.cpp -o $@

$(BUILD_DIR)/Directions.o: Directions.cpp Directions.hpp $(BUILD_DIR)/Direction.o $(BUILD_DIR)/Mesh.o $(BUILD_DIR)/Random_Pickup.o $(BUILD_DIR)/RNG.o
	$(COMPILE) Directions.cpp -o $@

$(BUILD_DIR)/Display.o: Display.cpp Display.hpp $(BUILD_DIR)/utils.o
	$(COMPILE) Display.cpp -o $@

$(BUILD_DIR)/Double.o: Double.cpp Double.hpp $(BUILD_DIR)/Exception.o $(BUILD_DIR)/Display.o
	$(COMPILE) Double.cpp -o $@

$(BUILD_DIR)/Eval_Point.o: Eval_Point.cpp Eval_Point.hpp  $(BUILD_DIR)/Parameters.o $(BUILD_DIR)/Cache_File_Point.o \
              Set_Element.hpp
	$(COMPILE) Eval_Point.cpp -o $@

$(BUILD_DIR)/Evaluator.o: Evaluator.cpp Evaluator.hpp $(BUILD_DIR)/Priority_Eval_Point.o $(BUILD_DIR)/Stats.o
	$(COMPILE) Evaluator.cpp -o $@

$(BUILD_DIR)/Evaluator_Control.o: Evaluator_Control.cpp Evaluator_Control.hpp \
            $(BUILD_DIR)/Barrier.o $(BUILD_DIR)/Pareto_Front.o $(BUILD_DIR)/Slave.o $(BUILD_DIR)/Quad_Model.o        
	$(COMPILE) Evaluator_Control.cpp -o $@

$(BUILD_DIR)/Exception.o: Exception.cpp Exception.hpp
	$(COMPILE) Exception.cpp -o $@

$(BUILD_DIR)/Extended_Poll.o: Extended_Poll.cpp Extended_Poll.hpp Signature_Element.hpp \
                 Set_Element.hpp $(BUILD_DIR)/Mads.o
	$(COMPILE) Extended_Poll.cpp -o $@

$(BUILD_DIR)/L_Curve.o: L_Curve.cpp L_Curve.hpp $(BUILD_DIR)/Double.o Uncopyable.hpp
	$(COMPILE) L_Curve.cpp -o $@

$(BUILD_DIR)/LH_Search.o: LH_Search.cpp LH_Search.hpp Search.hpp $(BUILD_DIR)/Mads.o $(BUILD_DIR)/RNG.o $(BUILD_DIR)/Evaluator_Control.o
	$(COMPILE) LH_Search.cpp -o $@

$(BUILD_DIR)/Mads.o: Mads.cpp Mads.hpp $(BUILD_DIR)/Evaluator_Control.o $(BUILD_DIR)/L_Curve.o \
        LH_Search.hpp LH_Search.cpp \
        Speculative_Search.cpp Speculative_Search.hpp \
        Extended_Poll.cpp Extended_Poll.hpp \
        VNS_Search.hpp VNS_Search.cpp \
        Quad_Model_Search.hpp Quad_Model_Search.cpp \
        Cache_Search.hpp Cache_Search.cpp \
        Phase_One_Search.cpp Phase_One_Search.hpp
	$(COMPILE) Mads.cpp -o $@

$(BUILD_DIR)/Mesh.o: Mesh.cpp Mesh.hpp $(BUILD_DIR)/Point.o
	$(COMPILE) Mesh.cpp -o $@

$(BUILD_DIR)/Multi_Obj_Evaluator.o: Multi_Obj_Evaluator.cpp Multi_Obj_Evaluator.hpp $(BUILD_DIR)/Phase_One_Evaluator.o
	$(COMPILE) Multi_Obj_Evaluator.cpp -o $@

$(BUILD_DIR)/Model_Sorted_Point.o: Model_Sorted_Point.cpp Model_Sorted_Point.hpp $(BUILD_DIR)/Point.o
	$(COMPILE) Model_Sorted_Point.cpp -o $@

$(BUILD_DIR)/Model_Stats.o: Model_Stats.cpp Model_Stats.hpp $(BUILD_DIR)/Double.o
	$(COMPILE) Model_Stats.cpp -o $@

$(BUILD_DIR)/nomad.o: nomad.cpp nomad.hpp  $(BUILD_DIR)/Mads.o 
	$(COMPILE) nomad.cpp -o $@

$(BUILD_DIR)/Parameters.o: Parameters.cpp Parameters.hpp  $(BUILD_DIR)/Parameter_Entries.o $(BUILD_DIR)/Signature.o
	$(COMPILE) Parameters.cpp -o $@

$(BUILD_DIR)/Parameter_Entries.o: Parameter_Entries.cpp Parameter_Entries.hpp $(BUILD_DIR)/Parameter_Entry.o
	$(COMPILE) Parameter_Entries.cpp -o $@

$(BUILD_DIR)/Parameter_Entry.o: Parameter_Entry.hpp Parameter_Entry.cpp  $(BUILD_DIR)/Display.o Uncopyable.hpp
	$(COMPILE) Parameter_Entry.cpp -o $@

$(BUILD_DIR)/Pareto_Front.o: Pareto_Front.cpp Pareto_Front.hpp $(BUILD_DIR)/Pareto_Point.o
	$(COMPILE) Pareto_Front.cpp -o $@

$(BUILD_DIR)/Pareto_Point.o: Pareto_Point.cpp Pareto_Point.hpp $(BUILD_DIR)/Multi_Obj_Evaluator.o
	$(COMPILE) Pareto_Point.cpp -o $@

$(BUILD_DIR)/Phase_One_Evaluator.o: Phase_One_Evaluator.cpp Phase_One_Evaluator.hpp $(BUILD_DIR)/Evaluator.o
	$(COMPILE) Phase_One_Evaluator.cpp -o $@

$(BUILD_DIR)/Phase_One_Search.o: Phase_One_Search.cpp Phase_One_Search.hpp $(BUILD_DIR)/Mads.o \
                    Search.hpp $(BUILD_DIR)/Evaluator_Control.o
	$(COMPILE) Phase_One_Search.cpp -o $@

$(BUILD_DIR)/Point.o: Point.cpp Point.hpp $(BUILD_DIR)/Double.o
	$(COMPILE) Point.cpp -o $@

$(BUILD_DIR)/Priority_Eval_Point.o: Priority_Eval_Point.cpp Priority_Eval_Point.hpp $(BUILD_DIR)/Eval_Point.o \
                       Set_Element.hpp
	$(COMPILE) Priority_Eval_Point.cpp -o $@

$(BUILD_DIR)/Quad_Model.o: Quad_Model.cpp Quad_Model.hpp $(BUILD_DIR)/Cache.o $(BUILD_DIR)/Model_Sorted_Point.o
	$(COMPILE) Quad_Model.cpp -o $@

$(BUILD_DIR)/Quad_Model_Evaluator.o: Quad_Model_Evaluator.cpp Quad_Model_Evaluator.hpp $(BUILD_DIR)/Evaluator.o Search.hpp
	$(COMPILE) Quad_Model_Evaluator.cpp -o $@

$(BUILD_DIR)/Quad_Model_Search.o: Quad_Model_Search.cpp Quad_Model_Search.hpp $(BUILD_DIR)/Mads.o $(BUILD_DIR)/Quad_Model_Evaluator.o
	$(COMPILE) Quad_Model_Search.cpp -o $@

$(BUILD_DIR)/Random_Pickup.o: Random_Pickup.cpp Random_Pickup.hpp RNG.cpp RNG.hpp Uncopyable.hpp
	$(COMPILE) Random_Pickup.cpp -o $@

$(BUILD_DIR)/RNG.o: RNG.cpp RNG.hpp defines.hpp
	$(COMPILE) RNG.cpp -o $@

$(BUILD_DIR)/Signature.o: Signature.cpp Signature.hpp $(BUILD_DIR)/Variable_Group.o
	$(COMPILE) Signature.cpp -o $@

$(BUILD_DIR)/Slave.o: Slave.cpp Slave.hpp $(BUILD_DIR)/Evaluator.o
	$(COMPILE) Slave.cpp -o $@

$(BUILD_DIR)/Speculative_Search.o: Speculative_Search.cpp Speculative_Search.hpp $(BUILD_DIR)/Mads.o Search.hpp $(BUILD_DIR)/Evaluator_Control.o
	$(COMPILE) Speculative_Search.cpp -o $@

$(BUILD_DIR)/Stats.o: Stats.cpp Stats.hpp $(BUILD_DIR)/Clock.o $(BUILD_DIR)/Double.o $(BUILD_DIR)/Model_Stats.o
	$(COMPILE) Stats.cpp -o $@

$(BUILD_DIR)/utils.o: utils.cpp utils.hpp defines.hpp 
	$(COMPILE) utils.cpp -o $@

$(BUILD_DIR)/Variable_Group.o: Variable_Group.cpp Variable_Group.hpp  $(BUILD_DIR)/Directions.o
	$(COMPILE) Variable_Group.cpp -o $@

$(BUILD_DIR)/VNS_Search.o: VNS_Search.cpp VNS_Search.hpp Search.hpp $(BUILD_DIR)/Evaluator_Control.o
	$(COMPILE) VNS_Search.cpp -o $@


