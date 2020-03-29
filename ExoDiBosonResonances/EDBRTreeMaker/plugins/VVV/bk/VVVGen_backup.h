#ifndef _VVVGen_backup_
#define _VVVGen_backup_
    // ************************* Gen Level Information******************//
    if(RunOnMC_)
    {//MC Info
        Int_t havegra=0;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {// loop on gen
            const reco::Candidate* ptop0 = &(*genParticles)[ik];
            const reco::Candidate* ptop=findLasttau(ptop0,6);
                if(ptop0->pdgId()== 6 && gentop_pt==-99) {
                    gentop_pt = ptop->pt();
                    gentop_eta = ptop->eta();
                    gentop_phi = ptop->phi();
                    gentop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        //if(abs(ptop->daughter(0)->pdgId())!=24&&abs(ptop->daughter(1)->pdgId())!=5) cout<<"no bW  "<<i<<"   "<<ptop->daughter(i)->pdgId()<<"   "<<ptop->daughter(i)->status()<<endl;
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            gent_w_pt=ptop->daughter(i)->pt();
                            gent_w_eta=ptop->daughter(i)->eta();
                            gent_w_phi=ptop->daughter(i)->phi();
                            gent_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                //if(abs(ptw->daughter(0)->pdgId())>=5) cout<<"no W-qq   "<<ptw->daughter(0)->pdgId()<<"   "<<ptw->daughter(1)->pdgId()<<endl;
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    gent_w_tag=4;
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) gent_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) gent_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) gent_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                        }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            gent_b_pt=ptop->daughter(i)->pt();
                            gent_b_eta=ptop->daughter(i)->eta();
                            gent_b_phi=ptop->daughter(i)->phi();
                            gent_b_mass=ptop->daughter(i)->mass();
                        }
                }
                }
                if(ptop0->pdgId()== -6 && genantitop_pt==-99) {
                    genantitop_pt = ptop->pt();
                    genantitop_eta = ptop->eta();
                    genantitop_phi = ptop->phi();
                    genantitop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        //cout<<i<<"   "<<ptop->daughter(i)->pdgId()<<"   "<<ptop->daughter(i)->status()<<endl;
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            genantit_w_pt=ptop->daughter(i)->pt();
                            genantit_w_eta=ptop->daughter(i)->eta();
                            genantit_w_phi=ptop->daughter(i)->phi();
                            genantit_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    genantit_w_tag=4;
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) genantit_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==12 ||abs(ptw->daughter(0)->pdgId())==13 ) genantit_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==15 ) genantit_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                            }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            genantit_b_pt=ptop->daughter(i)->pt();
                            genantit_b_eta=ptop->daughter(i)->eta();
                            genantit_b_phi=ptop->daughter(i)->phi();
                            genantit_b_mass=ptop->daughter(i)->mass();
                        }
                    }
                }
                
            //if((abs((*genParticles)[ik].pdgId())!=24)&&(abs((*genParticles)[ik].pdgId())!=9000024)&&(abs((*genParticles)[ik].pdgId())!=9000025))
            //cout<<"(*genParticles)[ik]->pdgId() "<<(*genParticles)[ik].pdgId()<<endl;
            if( abs((*genParticles)[ik].pdgId())==9000024|| abs((*genParticles)[ik].pdgId())==6) havegra=1;
        
            if( abs((*genParticles)[ik].pdgId())==9000024 || abs((*genParticles)[ik].pdgId())==6 )
            {//if Wkk
                const reco::Candidate* pwkk0 = &(*genParticles)[ik];
                const reco::Candidate* pwkk=findLastW(pwkk0,9000024);
                gen_gra_eta=pwkk->eta();
                gen_gra_m=pwkk->mass();
                gen_gra_pt=pwkk->pt();
                gen_gra_phi=pwkk->phi();
                for(int i=0;pwkk->daughter(i)!=NULL;i++)
                    {//loop on Wkk daughter
                       
                        if(abs(pwkk->daughter(i)->pdgId())==24){//if w
                            const reco::Candidate* pw0 = pwkk->daughter(i);
                            const reco::Candidate* pw= findLastW(pw0,24);                           //cout<<"check 4  "<<pw->daughter(0)->pdgId()<<endl;
                            if(pw->daughter(0)!=NULL)
                            {//loop on w daughter
                                const reco::Candidate* pl = pw->daughter(0);
                                if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                                {//beign of lep-w
                                    ptGenVlep = pw->pt();
                                    etaGenVlep = pw->eta();
                                    phiGenVlep = pw->phi();
                                    massGenVlep = pw->mass();
                                    status_1=0;

                                    for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                        const reco::Candidate* pl = pw->daughter(ii);
                                        if(abs(pl->pdgId())==11)
                                        {
                                            gen_ele_pt=pl->pt();
                                            gen_ele_eta=pl->eta();
                                            gen_ele_phi=pl->phi();
                                            gen_ele_e=pl->energy();
                                            status_1=1;
                                        }
                                        if(abs(pl->pdgId())==13)
                                        {
                                            gen_mu_pt=pl->pt();
                                            gen_mu_eta=pl->eta();
                                            gen_mu_phi=pl->phi();
                                            gen_mu_e=pl->energy();
                                            status_1=2;
                                        }
                                        if(abs(pl->pdgId())==15)
                                        {
                                            gen_tau_pt=pl->pt();
                                            gen_tau_eta=pl->eta();
                                            gen_tau_phi=pl->phi();
                                            gen_tau_e=pl->energy();
                                            const reco::Candidate* pl0= findLasttau(pl,15);
                                            for(int kk=0;pl0->daughter(kk)!=NULL&&kk<4;kk++)
                                            {
                                                pttau[kk]=pl0->daughter(kk)->pt();
                                                etatau[kk]=pl0->daughter(kk)->eta();
                                                phitau[kk]=pl0->daughter(kk)->pt();
                                                etau[kk]=pl0->daughter(kk)->energy();

                                            }
                                            status_1=3;
                                        }
                                        if(abs(pl->pdgId())==12)
                                        {
                                            gen_nele_pt=pl->pt();
                                            gen_nele_eta=pl->eta();
                                            gen_nele_phi=pl->phi();
                                            gen_nele_e=pl->energy();
                                        }
                                        if(abs(pl->pdgId())==14)
                                        {
                                            gen_nmu_pt=pl->pt();
                                            gen_nmu_eta=pl->eta();
                                            gen_nmu_phi=pl->phi();
                                            gen_nmu_e=pl->energy();
                                        }
                                        if(abs(pl->pdgId())==16)
                                        {
                                            gen_ntau_pt=pl->pt();
                                            gen_ntau_eta=pl->eta();
                                            gen_ntau_phi=pl->phi();
                                            gen_ntau_e=pl->energy();
                                        }

                                    }
                                }//end of if lep-w
                                numq=0;
                                for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                    const reco::Candidate* pl = pw->daughter(ii);
                                 
                                    if(abs(pl->pdgId())<6)
                                    {
                                        if(numq<3){
                                            ptq[numq]=pl->pt();
                                            etaq[numq]=pl->eta();
                                            phiq[numq]=pl->phi();
                                            eq[numq]=pl->energy();
                                            pdgidq[numq]=pl->pdgId();
                                            //cout<<numq<<"  "<<pdgidq[numq]<<endl;
                                            numq++;}
                                    }
                                }

                                if(abs(pl->pdgId())<6)
                                {
                                    ptGenVhad = pw->pt();
                                    etaGenVhad = pw->eta();
                                    phiGenVhad = pw->phi();
                                    massGenVhad = pw->mass();
                                    status_1=4;
                                }
                                if(abs(pl->pdgId())==24)
                                {
                                    status_1=5;
                                }
                                //if(status_1<0) cout<<"pw->daughter(0)  "<<pl->pdgId()<<endl;
                            }//end of loop on w daughter
                        }//end of if w
                    }//end of loop on Wkk daughter
            }//end of if Wkk

            //if(status_1<0) cout<<"(*genParticles)[ik].pdgId()"<<(*genParticles)[ik].pdgId()<<endl;
            if( havegra==0&&abs((*genParticles)[ik].pdgId())==24 )
            {//if W
                const reco::Candidate* pw0 = &(*genParticles)[ik];
                const reco::Candidate* pw=findFirstW(pw0,24);
                if(pw->mother(0)->pdgId()!=9000025){
                const reco::Candidate* pw1= findLastW(pw0,24);
                if(pw1->daughter(0)!=NULL){//loop on w daughter
                        const reco::Candidate* pl = pw1->daughter(0);

                    if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                    {//beign of lep-w
                        ptGenVlep = pw1->pt();
                        etaGenVlep = pw1->eta();
                        phiGenVlep = pw1->phi();
                        massGenVlep = pw1->mass();
                        status_1=0;

                        for(int ii=0;pw1->daughter(ii)!=NULL;ii++){
                            const reco::Candidate* pl = pw1->daughter(ii);
                            if(abs(pl->pdgId())==11)
                            {
                                gen_ele_pt=pl->pt();
                                gen_ele_eta=pl->eta();
                                gen_ele_phi=pl->phi();
                                gen_ele_e=pl->energy();
                                status_1=1;
                            }
                            if(abs(pl->pdgId())==13)
                            {
                                gen_mu_pt=pl->pt();
                                gen_mu_eta=pl->eta();
                                gen_mu_phi=pl->phi();
                                gen_mu_e=pl->energy();
                                status_1=2;
                            }
                            if(abs(pl->pdgId())==15)
                            {
                                gen_tau_pt=pl->pt();
                                gen_tau_eta=pl->eta();
                                gen_tau_phi=pl->phi();
                                gen_tau_e=pl->energy();
                                const reco::Candidate* pl0= findLasttau(pl,15);
                                for(int kk=0;pl0->daughter(kk)!=NULL&&kk<4;kk++)
                                {
                                    pttau[kk]=pl0->daughter(kk)->pt();
                                    etatau[kk]=pl0->daughter(kk)->eta();
                                    phitau[kk]=pl0->daughter(kk)->pt();
                                    etau[kk]=pl0->daughter(kk)->energy();
                                    pdgidtau[kk]=pl0->daughter(kk)->pdgId();

                                }
                                status_1=3;
                            }
                            if(abs(pl->pdgId())==12)
                            {
                                gen_nele_pt=pl->pt();
                                gen_nele_eta=pl->eta();
                                gen_nele_phi=pl->phi();
                                gen_nele_e=pl->energy();
                            }
                            if(abs(pl->pdgId())==14)
                            {
                                gen_nmu_pt=pl->pt();
                                gen_nmu_eta=pl->eta();
                                gen_nmu_phi=pl->phi();
                                gen_nmu_e=pl->energy();
                            }
                            if(abs(pl->pdgId())==16)
                            {
                                gen_ntau_pt=pl->pt();
                                gen_ntau_eta=pl->eta();
                                gen_ntau_phi=pl->phi();
                                gen_ntau_e=pl->energy();
                            }
                        }
                    
                    }//end of if lep-w
                        numq=0;
                        for(int ii=0;pw1->daughter(ii)!=NULL;ii++){
                            const reco::Candidate* pl = pw1->daughter(ii);
                            if(abs(pl->pdgId())<6)
                            {
                                if(numq<3){
                                    ptq[numq]=pl->pt();
                                    etaq[numq]=pl->eta();
                                    phiq[numq]=pl->phi();
                                    eq[numq]=pl->energy();
                                    pdgidq[numq]=pl->pdgId();
                                    numq++;}
                            }
                        }
                    if(abs(pl->pdgId())<6)
                    {
                        ptGenVhad = pw1->pt();
                        etaGenVhad = pw1->eta();
                        phiGenVhad = pw1->phi();
                        massGenVhad = pw1->mass();
                        status_1=4;

                    }
                    if(abs(pl->pdgId())==24)
                    {
                        status_1=5;
                    }
                }
            }
        }
        
		if( abs((*genParticles)[ik].pdgId())==9000025 ){//if Radion
            gen_rad_m=(*genParticles)[ik].mass();
            gen_rad_pt=(*genParticles)[ik].pt();
            gen_rad_phi=(*genParticles)[ik].phi();

            gen_rad_eta=(*genParticles)[ik].eta();
            for(int i=0;(*genParticles)[ik].daughter(i)!=NULL;i++){//loop on Radion daughter

                if(((*genParticles)[ik].daughter(i)->pdgId())==24){//if w-
                    const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                    const reco::Candidate* pw= findLastW(pw0,24);
                    if(pw->daughter(0)!=NULL){//loop on w daughter
                        const reco::Candidate* pl = pw->daughter(0);
                        if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)||(abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16)){//beign of lep-w
                            ptGenV_2 = pw->pt();
                            etaGenV_2 = pw->eta();
                            phiGenV_2 = pw->phi();
                            massGenV_2 = pw->mass();
                            status_2=0;

                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                                if(abs(pl->pdgId())==11)
                                {
                                    gen_ele_pt_2=pl->pt();
                                    gen_ele_eta_2=pl->eta();
                                    gen_ele_phi_2=pl->phi();
                                    gen_ele_e_2=pl->energy();
                                    status_2=1;
                                }
                                if(abs(pl->pdgId())==13)
                                {
                                    gen_mu_pt_2=pl->pt();
                                    gen_mu_eta_2=pl->eta();
                                    gen_mu_phi_2=pl->phi();
                                    gen_mu_e_2=pl->energy();
                                    status_2=2;
                                }
                                if(abs(pl->pdgId())==15)
                                {
                                    gen_tau_pt_2=pl->pt();
                                    gen_tau_eta_2=pl->eta();
                                    gen_tau_phi_2=pl->phi();
                                    gen_tau_e_2=pl->energy();
                                    const reco::Candidate* pl0= findLasttau(pl,15);
                                    for(int kk=0;pl0->daughter(kk)!=NULL&&kk<4;kk++)
                                    {
                                        pttau_2[kk]=pl0->daughter(kk)->pt();
                                        etatau_2[kk]=pl0->daughter(kk)->eta();
                                        phitau_2[kk]=pl0->daughter(kk)->pt();
                                        etau_2[kk]=pl0->daughter(kk)->energy();
                                        pdgidtau_2[kk]=pl0->daughter(kk)->pdgId();

                                    }

                                    status_2=3;
                                }
                                if(abs(pl->pdgId())==12)
                                {
                                    gen_nele_pt_2=pl->pt();
                                    gen_nele_eta_2=pl->eta();
                                    gen_nele_phi_2=pl->phi();
                                    gen_nele_e_2=pl->energy();
                                }
                                if(abs(pl->pdgId())==14)
                                {
                                    gen_nmu_pt_2=pl->pt();
                                    gen_nmu_eta_2=pl->eta();
                                    gen_nmu_phi_2=pl->phi();
                                    gen_nmu_e_2=pl->energy();
                                }
                                if(abs(pl->pdgId())==16)
                                {
                                    gen_ntau_pt_2=pl->pt();
                                    gen_ntau_eta_2=pl->eta();
                                    gen_ntau_phi_2=pl->phi();
                                    gen_ntau_e_2=pl->energy();
                                }

                            }
                        }//end of if lep-w
                             
                        numq_2=0;
                        for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                            const reco::Candidate* pl = pw->daughter(ii);
                            if(abs(pl->pdgId())<6)
                            {
                                if(numq_2<3){
                                    ptq_2[numq_2]=pl->pt();
                                    etaq_2[numq_2]=pl->eta();
                                    phiq_2[numq_2]=pl->phi();
                                    eq_2[numq_2]=pl->energy();
                                    pdgidq_2[numq_2]=pl->pdgId();
                                    numq_2++;}
                            }
                        }

                        if(abs(pl->pdgId())<6)
                            {
                                ptGenVhad_2 = pw->pt();
                                etaGenVhad_2 = pw->eta();
                                phiGenVhad_2 = pw->phi();
                                massGenVhad_2 = pw->mass();
                                status_2=4;
  
                            }
                        if(abs(pl->pdgId())==24)
                            {
                                status_2=5;
                            }
                    }//end of loop on w daughter
                }//end of if w-
                if(((*genParticles)[ik].daughter(i)->pdgId())==-24){//if w+
                    const reco::Candidate* pw0 = (*genParticles)[ik].daughter(i);
                    const reco::Candidate* pw= findLastW(pw0,24);
                    //cout<<((*genParticles)[ik].daughter(i)->pdgId())<<endl;
                        if(pw->daughter(0)!=NULL)
                        {//loop on w daughter
                            const reco::Candidate* pl = pw->daughter(0);
                            //cout<<(pl->pdgId())<<endl;
                            if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==13)|| (abs(pl->pdgId())==15)||(abs(pl->pdgId())==12) || (abs(pl->pdgId())==14)|| (abs(pl->pdgId())==16))
                            {//beign of lep-w
                                ptGenV_3 = pw->pt();
                                etaGenV_3 = pw->eta();
                                phiGenV_3 = pw->phi();
                                massGenV_3 = pw->mass();
                                status_3=0;
                            for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                                const reco::Candidate* pl = pw->daughter(ii);
                                if(abs(pl->pdgId())==11)
                                {
                                    gen_ele_pt_3=pl->pt();
                                    gen_ele_eta_3=pl->eta();
                                    gen_ele_phi_3=pl->phi();
                                    gen_ele_e_3=pl->energy();
                                    status_3=1;
                                }
                                if(abs(pl->pdgId())==13)
                                {
                                    gen_mu_pt_3=pl->pt();
                                    gen_mu_eta_3=pl->eta();
                                    gen_mu_phi_3=pl->phi();
                                    gen_mu_e_3=pl->energy();
                                    status_3=2;
                                }
                                if(abs(pl->pdgId())==15)
                                {
                                    gen_tau_pt_3=pl->pt();
                                    gen_tau_eta_3=pl->eta();
                                    gen_tau_phi_3=pl->phi();
                                    gen_tau_e_3=pl->energy();
                                    const reco::Candidate* pl0= findLasttau(pl,15);
                                    for(int kk=0;pl0->daughter(kk)!=NULL&&kk<4;kk++)
                                    {
                                        pttau_3[kk]=pl0->daughter(kk)->pt();
                                        etatau_3[kk]=pl0->daughter(kk)->eta();
                                        phitau_3[kk]=pl0->daughter(kk)->pt();
                                        etau_3[kk]=pl0->daughter(kk)->energy();
                                        pdgidtau_3[kk]=pl0->daughter(kk)->pdgId();
                                    }
                                    status_3=3;
                                }
                                if(abs(pl->pdgId())==12)
                                {
                                    gen_nele_pt_3=pl->pt();
                                    gen_nele_eta_3=pl->eta();
                                    gen_nele_phi_3=pl->phi();
                                    gen_nele_e_3=pl->energy();
                                }
                                if(abs(pl->pdgId())==14)
                                {
                                    gen_nmu_pt_3=pl->pt();
                                    gen_nmu_eta_3=pl->eta();
                                    gen_nmu_phi_3=pl->phi();
                                    gen_nmu_e_3=pl->energy();
                                }
                                if(abs(pl->pdgId())==16)
                                {
                                    gen_ntau_pt_3=pl->pt();
                                    gen_ntau_eta_3=pl->eta();
                                    gen_ntau_phi_3=pl->phi();
                                    gen_ntau_e_3=pl->energy();
                                }

                            }
                        }//end of if lep-w
                             
                        numq_3=0;
                        for(int ii=0;pw->daughter(ii)!=NULL;ii++){
                            const reco::Candidate* pl = pw->daughter(ii);
                            if(abs(pl->pdgId())<6)
                            {
                                if(numq_3<3){
                                    ptq_3[numq_3]=pl->pt();
                                    etaq_3[numq_3]=pl->eta();
                                    phiq_3[numq_3]=pl->phi();
                                    eq_3[numq_3]=pl->energy();
                                    pdgidq_3[numq_3]=pl->pdgId();
                                    numq_3++;}
                            }
                        }

                        if(abs(pl->pdgId())<6){
                            ptGenVhad_3 = pw->pt();
                            etaGenVhad_3 = pw->eta();
                            phiGenVhad_3 = pw->phi();
                            massGenVhad_3 = pw->mass();
                            status_3=4;
 
                        }
                        if(abs(pl->pdgId())==24){
                            status_3=5;
                        }
                    }//end of loop on w daughter
                }//end of if w+
			 
            }//end of loop on Radion daughter
 		}//end of if Radion

    }//end of loop on gen

    if(gen_mu_pt>0. && mus->size()>0 ){
        double drmumatch=10000.;  size_t mk=0;
        for(size_t ik=0; ik<mus->size();ik++)
        {
            double drtemp=deltaR(gen_mu_eta,gen_mu_phi,(*mus)[ik].eta(),(*mus)[ik].phi());
            if (drtemp<drmumatch) {drmumatch=drtemp; mk=ik;}
        }
        genmatch_mu_pt=(*mus)[mk].pt();
        genmatch_mu_eta=(*mus)[mk].eta();
        genmatch_mu_phi=(*mus)[mk].phi();
        genmatch_mu_e=(*mus)[mk].energy();
        genmatch_mu_dr=drmumatch;
    }
    if(gen_ele_pt>0. && eles->size()>0)
    {
        double drelematch=10000.;  size_t mk=0;
        for(size_t ik=0; ik<eles->size();ik++)
        {
            double drtemp=deltaR(gen_ele_eta,gen_ele_phi,(*eles)[ik].eta(),(*eles)[ik].phi());
            if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
        }
        genmatch_ele_pt=(*eles)[mk].pt();
        genmatch_ele_eta=(*eles)[mk].eta();
        genmatch_ele_phi=(*eles)[mk].phi();
        genmatch_ele_e=(*eles)[mk].energy();
        genmatch_ele_dr=drelematch;
        }
        
    //w and top info
        for( auto p=genParticles->begin(); p!= genParticles->end(); ++p)
        {}//std::cout<<p->pdgId()<<" "<<p->status()<<std::endl;}
        
        int igenw=0;
        int sizew=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==24)
                {
                    const reco::Candidate* pwtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pwtmp=findLastW(pwtmp1,24);
                    //const reco::Candidate* pwtmp2=findLasttau(pwtmp1,24);
                    //cout<<"genw"<<pwtmp->pt()<<"   "<<pwtmp1->pt()<<endl;
                    int woverlap=0;
                    for (int ia=0;ia<igenw;ia++){
                        if(pwtmp->pt()==ptgenwl[ia]) woverlap=1;
                    }
                    if(pwtmp->pt()>50&&igenw<sizew&&woverlap==0){
                    ptgenwl[igenw] = pwtmp->pt();
                    etagenwl[igenw] = pwtmp->eta();
                    phigenwl[igenw] = pwtmp->phi();
                    massgenwl[igenw] = pwtmp->mass();
                    const reco::Candidate* pwtmp2=findFirstW(pwtmp1,24);
                    ptgenwf[igenw] = pwtmp2->pt();
                    etagenwf[igenw] = pwtmp2->eta();
                    phigenwf[igenw] = pwtmp2->phi();
                    massgenwf[igenw] = pwtmp2->mass();
                        taggenwmother[igenw]=pwtmp2->mother(0)->pdgId();
                        /*cout<<ptgenwl[igenw]<<"   "<<ptgenwf[igenw]<<phigenwl[igenw]<<"   "<<phigenwf[igenw]<<"   "<<etagenwl[igenw]<<"   "<<etagenwf[igenw]<<" h0  "<<endl;
                        if(ptgenwl[igenw]!=ptgenwf[igenw] && pwtmp2->daughter(0)!=NULL) {cout<<phigenwl[igenw]<<"   "<<phigenwf[igenw]<<"   "<<etagenwl[igenw]<<"   "<<etagenwf[igenw]<<"  "<<pwtmp2->daughter(0)->pdgId()<<" h1  "<<endl;}
                        if(ptgenwl[igenw]!=ptgenwf[igenw] && pwtmp2->daughter(0)!=NULL&& pwtmp2->daughter(1)!=NULL) {cout<<phigenwl[igenw]<<"   "<<phigenwf[igenw]<<"   "<<etagenwl[igenw]<<"   "<<etagenwf[igenw]<<pwtmp2->daughter(0)->pdgId()<<"   ";
                            cout<<pwtmp2->daughter(1)->pdgId()<<endl;}*/
                    //for(int i=0;pw->daughter(i)!=NULL;i++)//loop on w daughter
                    if(pwtmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pwtmp->daughter(0);
                        //std::cout<< "pl pdgId" << pl->pdgId() << std::endl;
                         if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                                taggenwl[igenw]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenwl[igenw]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenwl[igenw]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenwl[igenw]=4;
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();
                            //cout<<"w_q  "<<igenw<<"  "<<pwtmp->daughter(0)->pt()<<"   "<<pwtmp->daughter(1)->pt()<<endl;
                            }
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ||(abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ||(abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16)){
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();

                        }
                    }
                    igenw+=1;
                    }
                }//end of if w

        }

        int igenz=0;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==23)
            {
                const reco::Candidate* pztmp1 = &(*genParticles)[ik];
                const reco::Candidate* pztmp=findLasttau(pztmp1,23);
                //const reco::Candidate* pztmp2=findLasttau(pztmp1,24);
                //cout<<pztmp->pt()<<"   "<<pztmp2->pt()<<endl;
                int zoverlap=0;
                for (int ia=0;ia<igenz;ia++){
                    if(pztmp->pt()==ptgenzl[ia]) zoverlap=1;}
                if(pztmp->pt()>50&&igenz<sizew&&zoverlap==0){
                    ptgenzl[igenz] = pztmp->pt();
                    etagenzl[igenz] = pztmp->eta();
                    phigenzl[igenz] = pztmp->phi();
                    massgenzl[igenz] = pztmp->mass();
                    const reco::Candidate* pztmp2=findFirstW(pztmp1,23);
                    ptgenzf[igenz] = pztmp2->pt();
                    etagenzf[igenz] = pztmp2->eta();
                    phigenzf[igenz] = pztmp2->phi();
                    massgenzf[igenz] = pztmp2->mass();
                    //for(int i=0;pz->daughter(i)!=NULL;i++)//loop on w daughter
                    if(pztmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pztmp->daughter(0);
                        //std::cout<< "pl pdgId" << pl->pdgId() << std::endl;
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                            taggenzl[igenz]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenzl[igenz]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenzl[igenz]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenzl[igenz]=4;}
                    }
                    igenz+=1;
                }
            }//end of if w
        }

        int igeng=0;
        int sizeg=10;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            //std::cout<<(*genParticles)[ik].pdgId()<<" "<<(*genParticles)[ik].status()<<std::endl;
            
            //		if( (*genParticles)[ik].pdgId()==5100039 ) // && (*genParticles)[ik].status()==3)//graviton
            //		{
            if(abs((*genParticles)[ik].pdgId())==21)
            {
                const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                const reco::Candidate* pgtmp=findLasttau(pgtmp1,21);
                int goverlap=0;
                for (int ia=0;ia<igeng;ia++){
                    if(pgtmp->pt()==ptgengl[ia]) goverlap=1;}
                if(pgtmp->pt()>50&&igeng<sizeg&&goverlap==0){
                    ptgengl[igeng] = pgtmp->pt();
                    etagengl[igeng] = pgtmp->eta();
                    phigengl[igeng] = pgtmp->phi();
                    egengl[igeng] = pgtmp->energy();
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp,21);
                    ptgengf[igeng] = pgtmp2->pt();
                    etagengf[igeng] = pgtmp2->eta();
                    phigengf[igeng] = pgtmp2->phi();
                    egengf[igeng] = pgtmp2->energy();
                    mothergengf[igeng] = pgtmp2->mother(0)->pdgId();
                    const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                    if (pgtmp3->mother(0)!=NULL) mmothergengf[igeng] = pgtmp3->mother(0)->pdgId();
                    //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergengf[igeng]<<endl;
                    igeng+=1;
                }
            }//end of if w
            
        }
     
        int igenq1=0;
        int sizeq1=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
                if(abs((*genParticles)[ik].pdgId())==1)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLasttau(pgtmp1,1);
                    int goverlap=0;
                    for (int ia=0;ia<igenq1;ia++){
                        if(pgtmp->pt()==ptgenq1l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp,1);
                    if(pgtmp->pt()>50&&igenq1<sizeq1&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq1l[igenq1] = pgtmp->pt();
                        etagenq1l[igenq1] = pgtmp->eta();
                        phigenq1l[igenq1] = pgtmp->phi();
                        egenq1l[igenq1] = pgtmp->energy();
                        ptgenq1f[igenq1] = pgtmp2->pt();
                        etagenq1f[igenq1] = pgtmp2->eta();
                        phigenq1f[igenq1] = pgtmp2->phi();
                        egenq1f[igenq1] = pgtmp2->energy();
                        mothergenq1f[igenq1] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq1f[igenq1] = pgtmp3->mother(0)->pdgId();
                        //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergenq1f[igenq1]<<endl;
                        igenq1+=1;
                        //cout<<"q1   "<<igenq1<<"   "<<pgtmp2->pdgId()<<"   "<<pgtmp->mother(0)->pdgId()<<"   "<<pgtmp->pt()<<"   "<<pgtmp2->pt()<<"   "<<genantit_w_q1_pt<<"  "<<genantit_w_q2_pt<<"   "<<gent_w_q1_pt<<"  "<<gent_w_q2_pt<<endl;
                    }
            }
        }

            int igenq2=0;
            int sizeq2=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==2)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLasttau(pgtmp1,2);
                    int goverlap=0;
                    for (int ia=0;ia<igenq2;ia++){
                        if(pgtmp->pt()==ptgenq2l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp1,2);
                    if(pgtmp->pt()>50&&igenq2<sizeq2&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq2l[igenq2] = pgtmp->pt();
                        etagenq2l[igenq2] = pgtmp->eta();
                        phigenq2l[igenq2] = pgtmp->phi();
                        egenq2l[igenq2] = pgtmp->energy();
                        ptgenq2f[igenq2] = pgtmp2->pt();
                        etagenq2f[igenq2] = pgtmp2->eta();
                        phigenq2f[igenq2] = pgtmp2->phi();
                        egenq2f[igenq2] = pgtmp2->energy();
                        mothergenq2f[igenq2] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq2f[igenq2] = pgtmp3->mother(0)->pdgId();
                        //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergenq2f[igenq2]<<endl;
                        igenq2+=1;
                    }
                }
            }

            int igenq3=0;
            int sizeq3=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==3)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLasttau(pgtmp1,3);
                    int goverlap=0;
                    for (int ia=0;ia<igenq3;ia++){
                        if(pgtmp->pt()==ptgenq3l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp1,3);
                    if(pgtmp->pt()>50&&igenq3<sizeq3&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq3l[igenq3] = pgtmp->pt();
                        etagenq3l[igenq3] = pgtmp->eta();
                        phigenq3l[igenq3] = pgtmp->phi();
                        egenq3l[igenq3] = pgtmp->energy();
                        ptgenq3f[igenq3] = pgtmp2->pt();
                        etagenq3f[igenq3] = pgtmp2->eta();
                        phigenq3f[igenq3] = pgtmp2->phi();
                        egenq3f[igenq3] = pgtmp2->energy();
                        mothergenq3f[igenq3] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq3f[igenq3] = pgtmp3->mother(0)->pdgId();
                        //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergenq3f[igenq3]<<endl;
                        igenq3+=1;
                    }
                }
            }

            int igenq4=0;
            int sizeq4=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==4)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pgtmp=findLasttau(pgtmp1,4);
                    int goverlap=0;
                    for (int ia=0;ia<igenq4;ia++){
                        if(pgtmp->pt()==ptgenq4l[ia]) goverlap=1;}
                    //if(pgtmp->pt()==gent_w_q1_pt||pgtmp->pt()==gent_w_q2_pt||pgtmp->pt()==genantit_w_q1_pt||pgtmp->pt()==genantit_w_q2_pt)
                    //    cout<<" q4 overlap with tw"<<endl;
                    //if(pgtmp->pt()==gent_b_pt||pgtmp->pt()==genantit_b_pt)                         cout<<" q4 overlap with tb"<<endl;
                    //if(pgtmp->pt()>50&&igenq4<sizeq4&&goverlap==0&&pgtmp->pt()!=gent_w_q1_pt&&pgtmp->pt()!=gent_w_q2_pt&&pgtmp->pt()!=genantit_w_q1_pt&&pgtmp->pt()!=genantit_w_q2_pt){
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp1,4);
                    if(pgtmp->pt()>50&&igenq4<sizeq4&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24){
                        ptgenq4l[igenq4] = pgtmp->pt();
                        etagenq4l[igenq4] = pgtmp->eta();
                        phigenq4l[igenq4] = pgtmp->phi();
                        egenq4l[igenq4] = pgtmp->energy();
                        ptgenq4f[igenq4] = pgtmp2->pt();
                        etagenq4f[igenq4] = pgtmp2->eta();
                        phigenq4f[igenq4] = pgtmp2->phi();
                        egenq4f[igenq4] = pgtmp2->energy();
                        mothergenq4f[igenq4] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq4f[igenq4] = pgtmp3->mother(0)->pdgId();
                        //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergenq4f[igenq4]<<endl;
                        igenq4+=1;
                    }
                }
            }

            int igenq5=0;
            int sizeq5=5;
            for(size_t ik=0; ik<genParticles->size();ik++)
            {
                if(abs((*genParticles)[ik].pdgId())==5)
                {
                    const reco::Candidate* pgtmp1 = &(*genParticles)[ik];
                    //const reco::Candidate* pgtmp=findLastW(pgtmp1,5);
                    const reco::Candidate* pgtmp=findLasttau(pgtmp1,5);
                    int goverlap=0;
                    for (int ia=0;ia<igenq5;ia++){
                        if(pgtmp->pt()==ptgenq5l[ia]) goverlap=1;}
                    const reco::Candidate* pgtmp2=findFirstW(pgtmp1,5);
                    //if(pgtmp->pt()==gent_w_q1_pt||pgtmp->pt()==gent_w_q2_pt||pgtmp->pt()==genantit_w_q1_pt||pgtmp->pt()==genantit_w_q2_pt)
                      //  cout<<" q5 overlap with tw"<<endl;
                    //if(pgtmp->pt()==gent_b_pt||pgtmp->pt()==genantit_b_pt)                         cout<<" q5 overlap with tb"<<endl;
                    if(pgtmp->pt()>50&&igenq5<sizeq5&&goverlap==0&&abs(pgtmp2->mother(0)->pdgId())!=24&&abs(pgtmp2->mother(0)->pdgId())!=6){
                        ptgenq5l[igenq5] = pgtmp->pt();
                        etagenq5l[igenq5] = pgtmp->eta();
                        phigenq5l[igenq5] = pgtmp->phi();
                        egenq5l[igenq5] = pgtmp->energy();
                        ptgenq5f[igenq5] = pgtmp2->pt();
                        etagenq5f[igenq5] = pgtmp2->eta();
                        phigenq5f[igenq5] = pgtmp2->phi();
                        egenq5f[igenq5] = pgtmp2->energy();
                        mothergenq5f[igenq5] = pgtmp2->mother(0)->pdgId();
                        const reco::Candidate* pgtmp3=findLasttau(pgtmp2->mother(0),pgtmp2->mother(0)->pdgId());
                        if (pgtmp3->mother(0)!=NULL) mmothergenq5f[igenq5] = pgtmp3->mother(0)->pdgId();
                        //cout<<pgtmp2->pdgId()<<"    "<<pgtmp2->mother(0)->pdgId()<<"  "<<mmothergenq5f[igenq5]<<endl;
                        igenq5+=1;
                        //cout<<"q5   "<<igenq5<<"   "<<pgtmp2->pdgId()<<"   "<<pgtmp->mother(0)->pdgId()<<"   "<<pgtmp->pt()<<"   "<<pgtmp2->pt()<<"   "<<gent_b_pt<<"   "<<genantit_b_pt<<endl;
                    }
                }
            }
        //cout<<"nng   "<<gentop_pt<<"  "<<genantitop_pt<<"  "<<igenw<<"  "<<igenq1<<"  "<<igenq2<<"  "<<igenq3<<"  "<<igenq4<<"  "<<igenq5<<"  "<<igeng<<"  "<<endl;

    }//end of MC Info
    // *************************End of Gen Level Information******************//
#endif
