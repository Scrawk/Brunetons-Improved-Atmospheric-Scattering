using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace BrunetonsImprovedAtmosphere
{
    public class RotateLight : MonoBehaviour
    {

        public float speed = 5.0f;

        private Vector3 lastMousePos;

        void Update()
        {
            Vector3 delta = (lastMousePos - Input.mousePosition) * Time.deltaTime * speed;

            if (Input.GetMouseButton(0) && Input.GetKey(KeyCode.LeftControl))
            {
                transform.Rotate(new Vector3(-delta.y, -delta.x, 0));
            }

            lastMousePos = Input.mousePosition;
        }
    }
}
