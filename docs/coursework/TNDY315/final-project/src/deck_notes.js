module.exports = {
  s1: `SPEAKER: DORON COHEN (about 30 seconds)

Hello. We are Doron Cohen, Manuel Blanco, and Amit Mekel. Our project is called SOLANGE. It is a research platform that asks one simple question: when do we really need a quantum computer to study a cancer mutation?

I will open with the problem and the technology. Manuel will explain what earlier projects teach us, our management approach, and the team. Amit will walk you through the plan, step by step. I will close.`,

  s2: `SPEAKER: DORON COHEN (about 1.5 minutes)

Lung cancer is the most common cancer in the world, and it kills more people than any other cancer. In 2022 there were almost 2.5 million new cases and about 1.8 million deaths.

Most patients have non-small cell lung cancer. Many of them carry mutations that no approved targeted drug can treat. My dissertation estimates that this applies to about 60 percent of these patients.

One reason is computational. To design a drug for a mutated protein, we need an accurate model of the electrons at the mutated site. When many electrons interact strongly, classical computers lose accuracy. This is where quantum computing enters the story.`,

  s3: `SPEAKER: DORON COHEN (about 1.5 minutes)

Quantum computers are expected to handle this kind of chemistry better. But today's machines are noisy. IBM plans its first large-scale, fault-tolerant machine, called Starling, for 2029. Our defense is planned for December 2027. So the technology we need matures after our milestones are due.

On top of that, this is what Rittel and Webber call a wicked problem. We cannot fully define what the platform must compute until we try to compute it. And every discipline on the team, chemistry, physics, biology, medicine, finance, and diplomacy, defines success in a different way.`,

  s4: `SPEAKER: DORON COHEN (about 1.5 minutes)

So SOLANGE does not assume that a quantum computer is needed. It checks.

For each mutation, taken from real patient sequencing data, the platform runs two independent classical methods, called DMRG and SHCI. If both methods agree that the problem can be solved classically, the verdict is Class B. If both agree that it cannot, the verdict is Class A, which means quantum-necessary.

A verification layer called LEON seals every result in nine provenance fields, P1 to P9, and checks it again before anyone may rely on it. This follows the audit-trail logic of the FDA rule for electronic records, known as Part 11.

So far, every real biological target we measured was classically tractable. We report that openly, because it is a real result.

This is not only a plan on paper.`,

  s4b: `SPEAKER: DORON COHEN (about 1 minute)

Here is the platform as it runs today. On the left is the classifier map. Each point is a real target. The horizontal line is the boundary between Class B and Class A. So far, every real biological target sits below the line, in Class B.

On the upper right are small test jobs that ran on real IBM Heron quantum computers in July and August. They proved that the pipeline works on real hardware.

On the lower right is the audit trail. Every run is sealed and can be verified again.

Now Manuel will explain what earlier projects teach us.`,

  s5: `SPEAKER: MANUEL BLANCO (about 2 minutes)

Thank you, Doron. I look at this field the way a banker and an investor would, so I start with money and trust.

The best-known failure is IBM Watson for Oncology. MD Anderson Cancer Center spent 62 million dollars on a project built on Watson, and then canceled it. There was very little published evidence that Watson helped patients. The promises came before the proof, and the computer scientists and the doctors never fully understood each other.

Quantum computing carries a similar risk. In 2023, IBM announced a result that classical computers supposedly could not match. Within months, other researchers matched it on classical computers, and some of the calculations even ran on a mobile phone.

On the positive side, real progress came slowly and steadily. A lung cancer target called KRAS was considered almost impossible to drug. The key discovery came in 2013, and the drug was approved in 2021. That is eight years of patient work. And Cleveland Clinic and IBM built a ten-year partnership that also trains people.

For an investor, the lesson is simple: let evidence set the pace, not marketing.`,

  s6: `SPEAKER: MANUEL BLANCO (about 1.5 minutes)

So how do we manage a project like this? We use a model from the PMBOK Guide called Cynefin.

Cynefin sorts problems by how well we understand cause and effect. In a complicated situation, the questions are known, so experts can analyze and plan ahead. In a complex situation, nobody knows the answer in advance, so the team tries something, measures the result, and adjusts, in short cycles.

SOLANGE has both kinds of work. The classification stage uses proven classical methods, so we plan it with fixed milestones. The quantum stage is complex, so we work in short cycles and update the plan every quarter.

Using different approaches for different parts of one project is what our course calls tailoring.`,

  s7: `SPEAKER: MANUEL BLANCO (about 1.5 minutes)

Cynefin tells us how to handle uncertainty, but it does not tell us how to talk to people. The PMBOK table does not link it to stakeholders. So we add a second model: cross-cultural communication.

This model says that the background of the sender shapes the message, and the background of the listener shapes how it is understood. A chemist, a doctor, an investor, and a health minister will read the same result in four different ways. Each of them defines success differently.

So we set three rules. Every result states clearly whether it is established or still experimental. Uncertain results are explained in a live meeting, not only in a written report. And the same result must mean the same thing to every reader.

From my side, this is how we protect the trust of funders.`,

  s8: `SPEAKER: MANUEL BLANCO (about 1 minute)

Here is how we divided the work. Each of us leads the area of his own professional expertise.

Doron leads the project and all the technical and scientific work.

I lead funding and investor relations. U.S. federal funding is not open to this project, so we focus on binational science grants and an academic collaboration with IBM Research Israel. I also own the budget and the reporting to funders.

Amit leads international partnerships and health diplomacy. He builds relationships with partner hospitals and health ministries, and he helps us find experts for the evaluation panel.`,

  s9: `SPEAKER: MANUEL BLANCO (about 1 minute)

This is our full timeline, from March 2026 to February 2028.

The upper part shows the research phases. The lower part shows the management tracks: the quarterly review of the plan, my funding track, and Amit's international partnership track.

There are four milestones, shown as diamonds: the proposal, the decision on access to IBM hardware, the freeze of our evaluation criteria, and the defense in December 2027. The red line shows where we are today.

Amit will now walk you through the plan, step by step.`,

  s10: `SPEAKER: AMIT MEKEL (about 1.5 minutes)

Thank you, Manuel. I will walk through the plan in simple terms.

Phases one and two are complete: the research foundation and a working prototype of the platform.

Phase 3A tests the design on a simulator. A simulator is a normal computer that imitates a quantum computer. The main test run finished on August 4, 2026. A larger run is still pending.

In July and August, the team also ran small test jobs on three real IBM quantum computers, called Kingston, Marrakesh, and Fez. These jobs proved that the whole pipeline works from end to end. They did not yet produce medical-grade chemistry results, and the team says that clearly. Being honest about what a result proves is part of the plan.`,

  s11: `SPEAKER: AMIT MEKEL (about 1.5 minutes)

Phase 3B is the hardest part of the plan. Here the project moves to IBM's newer machine, called Heron r3.

There are three problems. First, access depends on IBM, and it can take three to six months. Second, the hardware is still immature, so a run can fail for reasons the team does not control. Third, each discipline defines success differently, and here the differences are largest.

The plan answers each problem. Phase 3B is officially contingent. If access does not come by December, it becomes future work, and the project continues without delay. A target goes to the hardware only if it fits what the team has already proven it can run. And every discipline writes down its success criteria before the run, not after it.

In diplomacy we know this pattern well: when a result depends on another party, you prepare the alternative before the meeting, not after it.`,

  s12: `SPEAKER: AMIT MEKEL (about 1.5 minutes)

Phase 4, from January to April 2027, is the independent evaluation. It has six steps: a technical test system, comparative benchmarks, recruiting a panel of four to six experts, two rounds of panel sessions, walkthroughs of real cases, and a final report.

The two rounds matter. In the first round, the experts review the idea. In the second round, they use the platform themselves.

The evaluation criteria are frozen in advance on the Open Science Framework. This means nobody can change the rules after seeing the results.

The hardest step here is recruiting the experts. We need people from several fields, and we have only three weeks. This is where my international network helps, through partner hospitals and research institutions.`,

  s13: `SPEAKER: AMIT MEKEL (about 1.5 minutes)

Three management tracks run next to the research.

First, every quarter the whole team updates the plan based on what was actually achieved, never on vendor promises. This rule came from our earlier A3 analysis.

Second, Manuel's funding track runs from October 2026 to June 2027.

Third, my track builds partnerships with hospitals and health ministries. For me, this is health diplomacy. Cancer does not stop at borders. A shared and trusted scientific platform can bring countries together, even when politics divides them. After the results are final, we plan to share them with our partner health systems.`,

  s14: `SPEAKER: AMIT MEKEL (about 1.5 minutes)

Finally, the artifacts the project would use. These are the documents and tools that keep everyone aligned.

The project charter and the roadmap set the direction. The assumption log records every modeling choice behind a result. The risk register tracks the main risks: access, hardware, funding, and recruitment. The stakeholder register lists what success means for each group.

The change log records every formal change. For example, on September 23 the team retracted a result after a check found an error in the input structure. The change log shows that openly.

The communications plan and the budget complete the set.

Doron will now close.`,

  s15: `SPEAKER: DORON COHEN (about 1.5 minutes)

Thank you, Amit. Three lessons from this course shaped our plan.

First, tailoring. One project can use different approaches for different parts. Our classification stage follows a fixed plan, and our quantum stage works in short cycles.

Second, choose models for coverage. Cynefin handles the uncertainty of the technology, and the communication model handles the people. Each model covers what the other leaves open.

Third, evidence first. Our platform is allowed to say that a quantum computer is not needed. That is a result, not a failure.

With this structure, the uncertain future of quantum hardware becomes a managed condition. Every phase delivers value on its own, the hardware phase is contingent, and every result reaches its readers with its evidence and its context.`,

  s16: `SPEAKER: DORON COHEN (about 20 seconds)

Thank you for listening. On behalf of Manuel, Amit, and myself, we are happy to take your questions.`,
};
