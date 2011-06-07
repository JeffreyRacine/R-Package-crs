LIB                    = ../nomad.a
LIB_MPI                = ../nomad.MPI.a
COMPILATOR             = g++ -arch x86_64
COMPILATOR_MPI         = mpic++ -arch x86_64
COMPILATOR_OPTIONS     = -O3 -ansi -pedantic -DGCC_X
COMPILATOR_OPTIONS_MPI = -O3 -ansi -DGCC_X -DUSE_MPI
LIBS                   = -lm
INCLUDE                = -I.
COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
COMPILE_MPI            = $(COMPILATOR_MPI) $(COMPILATOR_OPTIONS_MPI) $(INCLUDE) -c
OBJS_LIB               = Barrier.o Cache.o Cache_File_Point.o Cache_Point.o \
                         Cache_Search.o Clock.o Direction.o Directions.o Display.o \
                         Double.o Eval_Point.o Evaluator.o Evaluator_Control.o \
                         Exception.o Extended_Poll.o L_Curve.o LH_Search.o Mads.o \
                         Mesh.o Model_Search.o Model_Stats.o Multi_Obj_Evaluator.o \
                         Parameters.o Parameter_Entries.o Parameter_Entry.o \
                         Pareto_Front.o Pareto_Point.o Phase_One_Evaluator.o \
                         Phase_One_Search.o Point.o Priority_Eval_Point.o Quad_Model.o \
                         Quad_Model_Evaluator.o Quad_Model_Sorted_Point.o \
                         Random_Pickup.o Signature.o Slave.o Speculative_Search.o \
                         Stats.o utils.o Variable_Group.o VNS_Search.o
OBJS                   = nomad.o $(OBJS_LIB)

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

all: clean $(OBJS)
	@echo "building the release version (no MPI)"
	@echo "exe file: "$(EXE)
	@$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)
	@strip $(EXE)
	@echo "building the library (no MPI)"
	@echo "lib file: "$(LIB)
	@ar -r $(LIB) $(OBJS_LIB)
	@ranlib $(LIB)

allmpi: COMPILE = $(COMPILE_MPI)
allmpi: clean $(OBJS)
	@echo "building the release version (MPI)"
	@echo "exe file: "$(EXE_MPI)
	@$(COMPILATOR_MPI) -o $(EXE_MPI) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS_MPI)
	@strip $(EXE_MPI)
	@echo "building the MPI library"
	@echo "lib file: "$(LIB_MPI)
	@ar -r $(LIB_MPI) $(OBJS_LIB)
	@ranlib $(LIB_MPI)

debug: COMPILATOR_OPTIONS = -ansi -Wall -DGCC_X -g -DDEBUG
debug: clean $(OBJS)
	@echo "building the debug version (no MPI)"
	@echo "exe file: "$(EXE)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

mpi: COMPILE = $(COMPILE_MPI)
mpi: clean $(OBJS)
	@echo "building the release version (MPI)"
	@echo "exe file: "$(EXE_MPI)
	@$(COMPILATOR_MPI) -o $(EXE_MPI) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS_MPI)
	@strip $(EXE_MPI)

lib: clean $(OBJS)
	@echo "building the library (no MPI)"
	@echo "lib file: "$(LIB)
	@ar -r $(LIB) $(OBJS_LIB)
	@ranlib $(LIB)

mpilib: COMPILE = $(COMPILE_MPI)
mpilib: clean $(OBJS)
	@echo "building the MPI library"
	@echo "lib file: "$(LIB_MPI)
	@ar -r $(LIB_MPI) $(OBJS_LIB)
	@ranlib $(LIB_MPI)

Barrier.o: Barrier.cpp Barrier.hpp Filter_Point.hpp Set_Element.hpp Eval_Point.o Cache.o
	$(COMPILE) Barrier.cpp

Cache.o: Cache.cpp Cache.hpp Cache_Point.o Clock.o
	$(COMPILE) Cache.cpp

Cache_File_Point.o: Cache_File_Point.cpp Cache_File_Point.hpp Display.o \
                    Eval_Point.cpp Eval_Point.hpp Uncopyable.hpp
	$(COMPILE) Cache_File_Point.cpp

Cache_Point.o: Cache_Point.cpp Cache_Point.hpp Eval_Point.o
	$(COMPILE) Cache_Point.cpp

Cache_Search.o: Cache_Search.cpp Cache_Search.hpp Mads.o Search.hpp Evaluator_Control.o
	$(COMPILE) Cache_Search.cpp

Clock.o: Clock.cpp Clock.hpp
	$(COMPILE) Clock.cpp

Direction.o: Direction.cpp Direction.hpp Point.o
	$(COMPILE) Direction.cpp

Directions.o: Directions.cpp Directions.hpp Direction.o Mesh.o Random_Pickup.o
	$(COMPILE) Directions.cpp

Display.o: Display.cpp Display.hpp utils.o
	$(COMPILE) Display.cpp

Double.o: Double.cpp Double.hpp Exception.o Display.o
	$(COMPILE) Double.cpp

Eval_Point.o: Eval_Point.cpp Eval_Point.hpp Parameters.o Cache_File_Point.o \
              Set_Element.hpp
	$(COMPILE) Eval_Point.cpp

Evaluator.o: Evaluator.cpp Evaluator.hpp Priority_Eval_Point.o Stats.o
	$(COMPILE) Evaluator.cpp

Evaluator_Control.o: Evaluator_Control.cpp Evaluator_Control.hpp \
                     Barrier.o Pareto_Front.o Slave.o Quad_Model.o
	$(COMPILE) Evaluator_Control.cpp

Exception.o: Exception.cpp Exception.hpp
	$(COMPILE) Exception.cpp

Extended_Poll.o: Extended_Poll.cpp Extended_Poll.hpp Mads.o Signature_Element.hpp \
                 Set_Element.hpp
	$(COMPILE) Extended_Poll.cpp

L_Curve.o: L_Curve.cpp L_Curve.hpp Double.o Uncopyable.hpp
	$(COMPILE) L_Curve.cpp

LH_Search.o: LH_Search.cpp LH_Search.hpp Mads.o Search.hpp Evaluator_Control.o
	$(COMPILE) LH_Search.cpp

Mads.o: Mads.cpp Mads.hpp Evaluator_Control.o L_Curve.o \
        LH_Search.hpp LH_Search.cpp \
        Speculative_Search.cpp Speculative_Search.hpp \
        Extended_Poll.cpp Extended_Poll.hpp \
        VNS_Search.hpp VNS_Search.cpp \
        Model_Search.hpp Model_Search.cpp \
        Cache_Search.hpp Cache_Search.cpp \
        Phase_One_Search.cpp Phase_One_Search.hpp
	$(COMPILE) Mads.cpp

Mesh.o: Mesh.cpp Mesh.hpp Point.o
	$(COMPILE) Mesh.cpp

Multi_Obj_Evaluator.o: Multi_Obj_Evaluator.cpp Multi_Obj_Evaluator.hpp \
                       Phase_One_Evaluator.o
	$(COMPILE) Multi_Obj_Evaluator.cpp

Model_Search.o: Model_Search.cpp Model_Search.hpp Mads.o Evaluator_Control.o \
                Quad_Model_Evaluator.o
	$(COMPILE) Model_Search.cpp

Model_Stats.o: Model_Stats.cpp Model_Stats.hpp Double.o
	$(COMPILE) Model_Stats.cpp

nomad.o: nomad.cpp nomad.hpp Model_Search.o VNS_Search.o Speculative_Search.o LH_Search.o
	$(COMPILE) nomad.cpp

Parameters.o: Parameters.cpp Parameters.hpp Parameter_Entries.o Signature.o
	$(COMPILE) Parameters.cpp

Parameter_Entries.o: Parameter_Entries.cpp Parameter_Entries.hpp Parameter_Entry.o
	$(COMPILE) Parameter_Entries.cpp

Parameter_Entry.o: Parameter_Entry.hpp Parameter_Entry.cpp Display.o Uncopyable.hpp
	$(COMPILE) Parameter_Entry.cpp

Pareto_Front.o: Pareto_Front.cpp Pareto_Front.hpp Pareto_Point.o
	$(COMPILE) Pareto_Front.cpp

Pareto_Point.o: Pareto_Point.cpp Pareto_Point.hpp Multi_Obj_Evaluator.o
	$(COMPILE) Pareto_Point.cpp

Phase_One_Evaluator.o: Phase_One_Evaluator.cpp Phase_One_Evaluator.hpp Evaluator.o
	$(COMPILE) Phase_One_Evaluator.cpp

Phase_One_Search.o: Phase_One_Search.cpp Phase_One_Search.hpp Mads.o \
                    Search.hpp Evaluator_Control.o
	$(COMPILE) Phase_One_Search.cpp

Point.o: Point.cpp Point.hpp Double.o
	$(COMPILE) Point.cpp

Priority_Eval_Point.o: Priority_Eval_Point.cpp Priority_Eval_Point.hpp Eval_Point.o \
                       Set_Element.hpp
	$(COMPILE) Priority_Eval_Point.cpp

Quad_Model.o: Quad_Model.cpp Quad_Model.hpp Cache.o Quad_Model_Sorted_Point.o
	$(COMPILE) Quad_Model.cpp

Quad_Model_Evaluator.o: Quad_Model_Evaluator.cpp Quad_Model_Evaluator.hpp \
                        Evaluator.o Search.hpp
	$(COMPILE) Quad_Model_Evaluator.cpp

Quad_Model_Sorted_Point.o: Quad_Model_Sorted_Point.cpp Quad_Model_Sorted_Point.hpp \
                           Eval_Point.o
	$(COMPILE) Quad_Model_Sorted_Point.cpp

Random_Pickup.o: Random_Pickup.cpp Random_Pickup.hpp Uncopyable.hpp
	$(COMPILE) Random_Pickup.cpp

Signature.o: Signature.cpp Signature.hpp Variable_Group.o
	$(COMPILE) Signature.cpp

Slave.o: Slave.cpp Slave.hpp Evaluator.o
	$(COMPILE) Slave.cpp

Speculative_Search.o: Speculative_Search.cpp Speculative_Search.hpp Mads.o Search.hpp \
                      Evaluator_Control.o
	$(COMPILE) Speculative_Search.cpp

Stats.o: Stats.cpp Stats.hpp Clock.o Double.o Model_Stats.o
	$(COMPILE) Stats.cpp

utils.o: utils.cpp utils.hpp defines.hpp 
	$(COMPILE) utils.cpp

Variable_Group.o: Variable_Group.cpp Variable_Group.hpp Directions.o
	$(COMPILE) Variable_Group.cpp

VNS_Search.o: VNS_Search.cpp VNS_Search.hpp Mads.o Search.hpp Evaluator_Control.o
	$(COMPILE) VNS_Search.cpp

clean:
	@echo "cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "cleaning trash files"
	@rm -f core *~
	@echo "cleaning obj files"
	@rm -f $(OBJS)
	@echo "cleaning exe file"
	@rm -f $(EXE) $(EXE_MPI)
	@echo "cleaning lib files"
	@rm -f $(LIB) $(LIB_MPI)
