from pipe import pipe

class myPipe(pipe):
    def __init__(self, input):
        pipe.__init__(self, input)
        self.group = 1

    def f1_a2b(self, runFlag):
        return None
    def f1_b2c(self, runFlag):
        codeList = ['ls', '-l', ['*py', '*pyc']]
        self.execM(codeList)
        return None
    def q2_c2d(self, runFlag):
        codeList = ['echo', self.input + ' ' + self.output]
        self.execM(codeList)
        return None
    def g3_d2E(self, runFlag):
        return None
    def g3_E2F(self, runFlag):
        # self.saveInput_('anoterInput.F')
        output = self.addExtraWord(self.getOutput_(), 'extra')
        self.saveOutput(output)
        return None
    def n4_F2g(self, runFlag):
        return None
    def g3_F2x(self, runFlag):
        return None
    def g3_d2e(self, runFlag):
        return None
    def q2_d2j(self, runFlag):
        return None
    def q2_c2r(self, runFlag):
        return None
    def f1_c2t(self, runFlag):
        return None


p = myPipe('input.a')
(p
    .branch(False)
    .run(p.f1_a2b, True) #1
    .run(p.f1_b2c, True) #2

        .branch(False and p.group == 1)
        .run(p.q2_c2d, True) #3
            
            .branch().run(p.g3_d2E, True) # 4
            .run(p.g3_E2F, True) # 5

                .branch().run(p.n4_F2g, True)
                .stop() # 6

            .run(p.g3_F2x, True) # 7 ## -> Uses continuation of latest 3's output
            .stop() 
            
            .branch().run(p.g3_d2e, True) # 8 # -> branch off from 2
            .stop()

        .run(p.q2_d2j, True) # 9 # Cont 2
        .stop()

        .branch().run(p.q2_c2r, True) # 10 # Branch off from 1
        .stop()

    .run(p.f1_c2t, True) # 11 # Con't from f2
    .stop()

    .branch()
    .run(p.f1_a2b, True) #1
    .stop()
)

print(p.i)