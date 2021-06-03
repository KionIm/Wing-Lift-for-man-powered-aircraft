function [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd)

%original foil 15�R�͌W��
x = linspace(2,6,41);
y = [200000,250000,300000,350000,400000,450000,500000,550000,600000,650000,700000];
z = zeros(11,41);
z(1,:) = [0.0174 0.01721 0.01703	0.01687	0.01673	0.01717	0.01745	0.01746	0.01739	0.01727	0.01712	0.01695	0.01679	0.01663	0.01647	0.01633	0.0162	0.01623	0.01657	0.01665	0.0166	0.01648	0.01631	0.01614	0.01596	0.0158	0.01567	0.01591	0.01602	0.01603	0.01597	0.01587	0.01575	0.01562	0.01565	0.01583	0.01591	0.01592	0.01587	0.01579	0.01572];
z(2,:) = [0.01338	0.01339	0.01352	0.01359	0.01361	0.01359	0.01354	0.01347	0.0134	0.01333	0.01326	0.0132	0.01315	0.01314	0.01326	0.01333	0.01335	0.01333	0.01327	0.0132	0.01313	0.01306	0.01302	0.01312	0.01319	0.01322	0.01323	0.01321	0.01317	0.01314	0.0132	0.0133	0.01336	0.0134	0.01341	0.01341	0.01344	0.01356	0.01366	0.01373	0.01377];
z(3,:) = [0.01152	0.01153	0.01151	0.01148	0.01144	0.0114	0.01137	0.01134	0.01132	0.01133	0.01138	0.01144	0.01147	0.01149	0.01148	0.01146	0.01143	0.0114	0.01138	0.01142	0.01148	0.01153	0.01155	0.01157	0.01157	0.01157	0.0116	0.01168	0.01174	0.0118	0.01183	0.01187	0.01191	0.01199	0.01208	0.01215	0.01222	0.01229	0.01237	0.01247	0.01257];
z(4,:) = [0.01015	0.01013	0.01011	0.01009	0.01008	0.01008	0.01012	0.01016	0.0102	0.01023	0.01025	0.01026	0.01026	0.01025	0.01024	0.01025	0.01029	0.01033	0.01037	0.0104	0.01042	0.01044	0.01047	0.01052	0.01058	0.01063	0.01068	0.01073	0.01078	0.01084	0.01091	0.01098	0.01105	0.01113	0.0112	0.01129	0.01139	0.01147	0.01156	0.01166	0.01177];
z(5,:) = [0.00917	0.00917	0.00918	0.00923	0.00926	0.00929	0.00931	0.00934	0.00935	0.00936	0.00937	0.00938	0.00939	0.00944	0.00947	0.0095	0.00954	0.00956	0.00959	0.00963	0.00967	0.00973	0.00978	0.00983	0.00988	0.00993	0.00999	0.01007	0.01013	0.0102	0.01027	0.01034	0.01043	0.01054	0.01061	0.01069	0.01078	0.01089	0.01101	0.01112	0.0112];
z(6,:) = [0.00852	0.00856	0.00858	0.0086	0.00863	0.00865	0.00867	0.00868	0.0087	0.00872	0.00875	0.00879	0.00882	0.00885	0.00888	0.00891	0.00895	0.00899	0.00905	0.0091	0.00914	0.00919	0.00924	0.0093	0.00937	0.00945	0.0095	0.00956	0.00963	0.00971	0.00981	0.0099	0.00997	0.01005	0.01014	0.01024	0.01037	0.01047	0.01055	0.01064	0.01075];
z(7,:) = [0.00803	0.00805	0.00808	0.0081	0.00811	0.00813	0.00816	0.00818	0.00823	0.00827	0.00829	0.00832	0.00835	0.00838	0.00842	0.00847	0.00854	0.00857	0.00862	0.00866	0.00871	0.00877	0.00885	0.00892	0.00897	0.00903	0.0091	0.00918	0.00927	0.00937	0.00943	0.00951	0.00959	0.0097	0.00982	0.00992	0.00999	0.01007	0.01018	0.01031	0.01045];
z(8,:) = [0.00762	0.00764	0.00766	0.00768	0.00771	0.00774	0.00779	0.00783	0.00785	0.00788	0.0079	0.00794	0.00798	0.00803	0.0081	0.00813	0.00817	0.00821	0.00826	0.00832	0.0084	0.00848	0.00852	0.00858	0.00864	0.00872	0.00881	0.00891	0.00897	0.00904	0.00912	0.00921	0.00933	0.00945	0.00951	0.00958	0.00968	0.0098	0.00995	0.01003	0.01012];
z(9,:) = [0.00728	0.0073	0.00733	0.00736	0.00741	0.00745	0.00747	0.0075	0.00752	0.00756	0.00759	0.00764	0.00771	0.00775	0.00779	0.00783	0.00787	0.00793	0.008	0.00809	0.00813	0.00818	0.00824	0.00831	0.00839	0.0085	0.00856	0.00863	0.0087	0.00879	0.00889	0.00902	0.00909	0.00915	0.00923	0.00934	0.00947	0.00959	0.00967	0.00977	0.00992];
z(10,:) = [0.00699	0.00703	0.00707	0.00713	0.00714	0.00717	0.00719	0.00722	0.00726	0.0073	0.00736	0.00742	0.00745	0.00749	0.00753	0.00758	0.00764	0.00772	0.00779	0.00783	0.00788	0.00795	0.00802	0.00812	0.00821	0.00826	0.00833	0.00841	0.0085	0.00862	0.00872	0.00877	0.00884	0.00893	0.00905	0.00919	0.00926	0.00935	0.00947	0.00963	0.00975];
z(11,:) = [0.00677	0.00682	0.00686	0.00688	0.0069	0.00693	0.00696	0.007	0.00705	0.00712	0.00716	0.00719	0.00722	0.00727	0.00733	0.0074	0.00748	0.00752	0.00757	0.00762	0.00769	0.00777	0.00788	0.00794	0.00799	0.00806	0.00815	0.00825	0.00837	0.00843	0.00849	0.00856	0.00866	0.00878	0.00891	0.00897	0.00907	0.0092	0.00937	0.00945	0.00956];

drag15 = interp2(x,y,z,alpha,Re,"spline");


drag15m = interp2(x,y,z,alpham,codemax*U/nu,"spline");


D = zeros(1,npartition);

for i = 1:npartition 
    
        D(1,i) = 1/2*drag15(1,i)*chord(1,i)*rho*U^2*nd;
        
end    






