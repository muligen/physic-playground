using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicShapes
{
    public class ThreeDBall
    {
        public Vector3 vel;
        public Vector3 pos;
        public float radius;
        PhysicScene s;

        public ThreeDBall(Vector3 vel, Vector3 pos, float radius)
        {
            this.vel = vel;
            this.pos = pos;
            this.radius = radius;
            s = PhysicScene.Instance;
        }

        public void simulate(float timestamp)
        {
            vel += timestamp * new Vector3(0, s.gravity, 0);
            pos += timestamp * vel;

            //¼ÆËãÇ½±ÚÅö×²
            if (pos.y < s.plane_lower_bound + radius)
            {
                pos.y = s.plane_lower_bound + radius; vel.y = -vel.y;
            }
            if (pos.x < s.x_lower_bound + radius)
            {
                pos.x = s.x_lower_bound + radius; vel.x = -vel.x;
            }
            if (pos.x > s.x_upper_bound - radius)
            {
                pos.x = s.x_upper_bound - radius; vel.x = -vel.x;
            }
            if (pos.z < s.z_lower_bound + radius)
            {
                pos.z = s.z_lower_bound + radius; vel.z = -vel.z;
            }
            if (pos.z > s.z_upper_bound - radius)
            {
                pos.z = s.z_upper_bound - radius; vel.z = -vel.z;
            }
        }
    }
}